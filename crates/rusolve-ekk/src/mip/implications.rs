//! HiGHS source: mip/HighsImplications.h, mip/HighsImplications.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// C++ `static constexpr bool kSkipBadVbds = true;` (HighsImplications.cpp:157)
const K_SKIP_BAD_VBDS: bool = true;
/// C++ `static constexpr bool kUseDualsForBreakingTies = true;` (HighsImplications.cpp:158)
const K_USE_DUALS_FOR_BREAKING_TIES: bool = true;
/// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications::Implication (lines 26-29)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HighsImplications_Implication {
    pub lb: f64,
    pub ub: f64,
}

impl Default for HighsImplications_Implication {
    fn default() -> Self {
        Self {
            lb: -kHighsInf,
            ub: kHighsInf,
        }
    }
}

/// HiGHS: HighsImplications.h::HighsImplications (lines 23-195)
/// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications (lines 23-160)
///
/// Removed C++ members:
/// * C++: `const HighsMipSolver & mipsolver`: removed; methods that use it take `ms: &mut HighsMipSolver`.
#[derive(Clone, Default, Debug)]
pub struct HighsImplications {
    /// C++: `HighsInt nextCleanupCall`
    pub next_cleanup_call: HighsInt,
    /// C++: `std::vector<Implics> implications`
    pub implications: Vec<HighsImplications_Implics>,
    /// C++: `int64_t numImplications`
    pub num_implications: i64,
    /// C++: `int64_t numVarBounds`
    pub num_var_bounds: i64,
    /// C++: `int64_t maxVarBounds`
    pub max_var_bounds: i64,
    /// C++: `std::vector<HighsHashTree<HighsInt, VarBound>> vubs`
    pub vubs: Vec<HighsHashTree<HighsInt, HighsImplications_VarBound>>,
    /// C++: `std::vector<HighsHashTree<HighsInt, VarBound>> vlbs`
    pub vlbs: Vec<HighsHashTree<HighsInt, HighsImplications_VarBound>>,
    /// C++: `std::vector<HighsHashTree<HighsInt, HighsInt>> rowToVarBounds`
    pub row_to_var_bounds: Vec<HighsHashTree<HighsInt, HighsInt>>,
    /// C++: `std::vector<HighsSubstitution> substitutions`
    pub substitutions: Vec<HighsSubstitution>,
    /// C++: `std::vector<uint8_t> colsubstituted`
    pub colsubstituted: Vec<u8>,
    /// C++ `std::function<void(HighsInt, HighsInt, HighsInt, double)> storeLiftingOpportunity`, set only by `HPresolve::runProbing` (HPresolve.cpp:1685) to push into `liftingOpportunities`: `Some(buffer)` = set; the call appends `(col, implcol, ..., val)` and HPresolve drains the buffer after each `implications` call.
    pub store_lifting_opportunity: Option<Vec<(HighsInt, HighsInt, HighsInt, f64)>>,
    /// C++: `std::vector<HighsDomainChange> implicationsDown` (HiGHS #3244, mip_dual_fix_probing;
    /// written only while the global domain's dual-fix probing is active, else empty)
    pub implications_down: Vec<HighsDomainChange>,
    /// C++: `std::vector<HighsDomainChange> implicationsUp` (HiGHS #3244, mip_dual_fix_probing;
    /// written only while the global domain's dual-fix probing is active, else empty)
    pub implications_up: Vec<HighsDomainChange>,
    /// C++: `std::vector<HighsInt> dualFixProbingBinInds_` (HiGHS #3244, mip_dual_fix_probing)
    pub dual_fix_probing_bin_inds_: Vec<HighsInt>,
    /// C++: `std::vector<TentativeFixing> dualFixProbingBinFlags_` (HiGHS #3244,
    /// mip_dual_fix_probing; sized only when the option is on)
    pub dual_fix_probing_bin_flags_: Vec<HighsImplications_TentativeFixing>,
    /// C++: `std::vector<HighsHashTree<HighsInt, Implication>> implications` (HiGHS #3209,
    /// mip_keep_implications; sized 2*ncols only when keep on)
    pub implication_tree: Vec<HighsHashTree<HighsInt, HighsImplications_Implication>>,
    /// C++: `std::vector<HighsHashTree<HighsInt, bool>> reverseImplications` (HiGHS #3209,
    /// mip_keep_implications; sized ncols only when keep on)
    pub reverse_implications: Vec<HighsHashTree<HighsInt, bool>>,
    /// C++: `std::vector<uint8_t> hasProbed` (HiGHS #3209, mip_keep_implications;
    /// sized 2*ncols only when keep on)
    pub has_probed: Vec<bool>,
    /// HiGHS #3209: set at construction from ms.options_mip_.mip_keep_implications
    pub keep_x1b: bool,
    /// Port-only scratch for the `cleanup_varbounds` snapshots ([`collect_var_bounds_into`]).
    /// Not part of the C++ state; always left empty or stale, never read across calls.
    pub vb_snapshot_: Vec<(HighsInt, HighsImplications_VarBound)>,
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::TentativeFixing::Direction (line 76)
/// C++: `enum Direction : uint8_t { Undecided, FixLower, FixUpper }`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HighsImplications_TentativeFixingDirection {
    #[default]
    Undecided,
    FixLower,
    FixUpper,
}

/// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::TentativeFixing (lines 75-94)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HighsImplications_TentativeFixing {
    /// C++: `Direction downProbe = Undecided`
    pub down_probe: HighsImplications_TentativeFixingDirection,
    /// C++: `Direction upProbe = Undecided`
    pub up_probe: HighsImplications_TentativeFixingDirection,
}

impl HighsImplications_TentativeFixing {
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::TentativeFixing::isUndecided (lines 80-82)
    pub fn is_undecided(&self) -> bool {
        self.down_probe == HighsImplications_TentativeFixingDirection::Undecided
            && self.up_probe == HighsImplications_TentativeFixingDirection::Undecided
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::TentativeFixing::record (lines 84-88)
    pub fn record(&mut self, up_probing: bool, boundtype: HighsBoundType) {
        let probe = if up_probing {
            &mut self.up_probe
        } else {
            &mut self.down_probe
        };
        if *probe != HighsImplications_TentativeFixingDirection::Undecided {
            return;
        }
        *probe = if boundtype == HighsBoundType::kLower {
            HighsImplications_TentativeFixingDirection::FixUpper
        } else {
            HighsImplications_TentativeFixingDirection::FixLower
        };
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::TentativeFixing::clear (lines 90-93)
    pub fn clear(&mut self) {
        self.down_probe = HighsImplications_TentativeFixingDirection::Undecided;
        self.up_probe = HighsImplications_TentativeFixingDirection::Undecided;
    }
}

/// HiGHS: HighsImplications.h::HighsImplications::Implics (lines 26-29)
#[derive(Clone, Debug)]
pub struct HighsImplications_Implics {
    /// C++: `std::vector<HighsDomainChange> implics`
    pub implics: Vec<HighsDomainChange>,
    /// C++: `bool computed = false`
    pub computed: bool,
}

impl Default for HighsImplications_Implics {
    fn default() -> Self {
        Self {
            implics: Default::default(),
            computed: false,
        }
    }
}

/// HiGHS: HighsImplications.h::HighsImplications::VarBound (lines 38-50)
/// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::VarBound (lines 51-68)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HighsImplications_VarBound {
    /// C++: `double coef`
    pub coef: f64,
    /// C++: `double constant`
    pub constant: f64,
    /// C++: `HighsInt origin`
    pub origin: HighsInt,
}

impl Default for HighsImplications_VarBound {
    fn default() -> Self {
        Self {
            coef: 0.0,
            constant: 0.0,
            origin: -1,
        }
    }
}

/// HiGHS: HighsImplications.h::HighsImplications::VarBound::minValue/maxValue (lines 42-49)
///
/// PORT NOTE: the C++ sums in `HighsCDouble` and rounds once to `double`; `HighsCDouble`
/// arithmetic is reproduced by `HighsCDouble::from(constant) + cpp_min(coef, 0.0)`.
fn var_bound_min_value(vb: &HighsImplications_VarBound) -> f64 {
    f64::from(HighsCDouble::from(vb.constant) + cpp_min(vb.coef, 0.0))
}

fn var_bound_max_value(vb: &HighsImplications_VarBound) -> f64 {
    f64::from(HighsCDouble::from(vb.constant) + cpp_max(vb.coef, 0.0))
}

/// HiGHS #3264 (mip_aggregate_implications): HighsImplications.h lines 51-68 at f2ae2b320 make
/// `VarBound::minValue/maxValue` return the plain sum for an infinite term instead of the
/// `HighsCDouble` NaN; every MIP user of the header (getBestVub/Vlb, addVUB/VLB,
/// cleanupVub/Vlb, HighsTransformedLp::transform) sees that in 1.16.0.
fn var_bound_min_value_x(vb: &HighsImplications_VarBound, inf_safe: bool) -> f64 {
    if inf_safe {
        var_bound_min_value_inf_safe(vb)
    } else {
        var_bound_min_value(vb)
    }
}

fn var_bound_max_value_x(vb: &HighsImplications_VarBound, inf_safe: bool) -> f64 {
    if inf_safe {
        var_bound_max_value_inf_safe(vb)
    } else {
        var_bound_max_value(vb)
    }
}

/// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::VarBound::minValue (lines 56-61)
fn var_bound_min_value_inf_safe(vb: &HighsImplications_VarBound) -> f64 {
    let m = cpp_min(vb.coef, 0.0);
    if vb.constant.abs() >= kHighsInf || m.abs() >= kHighsInf {
        return vb.constant + m;
    }
    f64::from(HighsCDouble::from(vb.constant) + m)
}

/// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::VarBound::maxValue (lines 62-67)
fn var_bound_max_value_inf_safe(vb: &HighsImplications_VarBound) -> f64 {
    let m = cpp_max(vb.coef, 0.0);
    if vb.constant.abs() >= kHighsInf || m.abs() >= kHighsInf {
        return vb.constant + m;
    }
    f64::from(HighsCDouble::from(vb.constant) + m)
}

/// Snapshot of one variable-bound tree as `(key, value)` pairs, in tree iteration order,
/// into `out` (cleared first).
///
/// PORT NOTE: `cleanupVarbounds` mutates each `VarBound` in its `for_each` lambda through
/// `cleanupVub`/`cleanupVlb`, which take `&mut self` here (they lend `self` back for the bound
/// change), so that lambda is written as snapshot, mutate, then `find_mut` write-back.
/// PORT NOTE (perf): the snapshot goes into the reused
/// `HighsImplications::vb_snapshot_` buffer instead of a fresh `Vec` per call.
fn collect_var_bounds_into(
    tree: &HighsHashTree<HighsInt, HighsImplications_VarBound>,
    out: &mut Vec<(HighsInt, HighsImplications_VarBound)>,
) {
    out.clear();
    tree.for_each(&mut |k, v| out.push((*k, *v)));
}

/// The `cleanupVub` (`is_vub`) / `cleanupVlb` lambda body of `cleanupVarbounds`
/// (HighsImplications.cpp:836-875 / 795-834, `allowBoundChanges = true`) run in place over
/// `tree` with `for_each_mut`, for every entry whose cleanup does not change a bound. Same
/// arithmetic, order and debug-solution checks as [`HighsImplications::cleanup_vub`] /
/// [`HighsImplications::cleanup_vlb`]; redundant keys are pushed to `del_vbds`.
///
/// Returns `None` when the whole tree was processed, or `Some(k)` when the entry at iteration
/// position `k` needs the bound change: that entry is left untouched (the bound-change branch
/// does not modify it) and nothing after it is visited, so the caller resumes at `k`.
fn cleanup_vbds_in_place(
    tree: &mut HighsHashTree<HighsInt, HighsImplications_VarBound>,
    ms: &HighsMipSolver,
    col: HighsInt,
    bound: f64,
    is_vub: bool,
    del_vbds: &mut Vec<HighsInt>,
) -> Option<usize> {
    let feastol = ms.mipdata_.feastol;
    let epsilon = ms.mipdata_.epsilon;
    let inf_safe = ms.options_mip_.mip_aggregate_implications;
    let debug_solution = &ms.mipdata_.debug_solution;
    let mut pos: usize = 0;
    let mut stop: Option<usize> = None;
    tree.for_each_mut(&mut |&vb_col, vb| {
        if stop.is_some() {
            return;
        }
        if vb_col == -1 {
            pos += 1;
            return;
        }
        let maxv = HighsCDouble::from(var_bound_max_value_x(vb, inf_safe));
        let minv = HighsCDouble::from(var_bound_min_value_x(vb, inf_safe));
        if is_vub {
            if minv >= bound - feastol {
                debug_solution.check_vub(col, vb_col, vb.coef, vb.constant);
                del_vbds.push(vb_col);
            } else if maxv > bound + epsilon {
                debug_solution.check_vub(col, vb_col, vb.coef, vb.constant);
                let newcoef = f64::from(HighsCDouble::from(bound) - minv);
                if vb.coef > 0.0 {
                    vb.coef = newcoef;
                } else {
                    vb.constant = bound;
                    vb.coef = -newcoef;
                }
                debug_solution.check_vub(col, vb_col, vb.coef, vb.constant);
            } else if maxv < bound - epsilon {
                stop = Some(pos);
                return;
            } else {
                debug_solution.check_vub(col, vb_col, vb.coef, vb.constant);
            }
        } else if maxv <= bound + feastol {
            debug_solution.check_vlb(col, vb_col, vb.coef, vb.constant);
            del_vbds.push(vb_col);
        } else if minv < bound - epsilon {
            debug_solution.check_vlb(col, vb_col, vb.coef, vb.constant);
            let newcoef = f64::from(HighsCDouble::from(bound) - maxv);
            if vb.coef < 0.0 {
                vb.coef = newcoef;
            } else {
                vb.constant = bound;
                vb.coef = -newcoef;
            }
            debug_solution.check_vlb(col, vb_col, vb.coef, vb.constant);
        } else if minv > bound + epsilon {
            stop = Some(pos);
            return;
        } else {
            debug_solution.check_vlb(col, vb_col, vb.coef, vb.constant);
        }
        pos += 1;
    });
    stop
}

/// Write a mutated `VarBound` back into its tree (see [`collect_var_bounds_into`]).
fn write_back_var_bound(
    tree: &mut HighsHashTree<HighsInt, HighsImplications_VarBound>,
    key: HighsInt,
    vbd: HighsImplications_VarBound,
) {
    if let Some(slot) = tree.find_mut(&key) {
        *slot = vbd;
    }
}

impl HighsImplications {
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::addImplication (lines 376-393)
    /// C++: `void addImplication(ImplIdx idx, HighsInt implCol, Implication implic)`
    pub fn add_implication(
        &mut self,
        col: HighsInt,
        val: bool,
        impl_col: HighsInt,
        implic: HighsImplications_Implication,
    ) {
        let idx = (2 * col + val as HighsInt) as usize;
        let (slot, inserted) =
            self.implication_tree[idx].insert_or_get(HighsHashTableEntry::new(impl_col, implic));
        if !inserted {
            if slot.lb == -kHighsInf && implic.lb != -kHighsInf {
                self.num_implications += 1;
            }
            if slot.ub == kHighsInf && implic.ub != kHighsInf {
                self.num_implications += 1;
            }
            slot.lb = cpp_max(implic.lb, slot.lb);
            slot.ub = cpp_min(implic.ub, slot.ub);
        } else {
            if implic.lb != -kHighsInf {
                self.num_implications += 1;
            }
            if implic.ub != kHighsInf {
                self.num_implications += 1;
            }
        }
        self.reverse_implications[impl_col as usize]
            .insert_or_get(HighsHashTableEntry::new(col, false));
    }

    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications::probedBefore (lines 120-122)
    /// C++: `bool probedBefore(const HighsInt col, const bool val) const`
    pub fn probed_before(&self, col: HighsInt, val: bool) -> bool {
        self.has_probed[(2 * col + val as HighsInt) as usize]
    }
    /// HiGHS: HighsImplications.cpp::HighsImplications::computeImplications (lines 15-152)
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::computeImplications (lines 15-211)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::computeImplications (lines 15-157)
    /// C++: `bool computeImplications(HighsInt col, bool val)`
    ///
    /// PORT NOTE: the C++ lambdas `storeLiftingOpportunities`, `doBacktrack` and `isInfeasible`
    /// are private helpers below (`store_lifting_opportunities`, `do_backtrack`,
    /// `is_infeasible_after_change`).
    ///
    /// PORT NOTE (HiGHS #3244): with the global domain's dual-fix probing inactive, the 1.15.1
    /// body runs unchanged. When active, the latest tail (lines 88-209) runs instead, with
    /// the 1.15.1 storage kept (the parts of #3244 not ported): `hasProbed[idx] = true`
    /// and the `addImplication` loop over the safe prefix become `implications[loc].computed
    /// = true` and the safe non-binary prefix stored in `implications[loc].implics`, sorted
    /// as 1.15.1 stores it (`pdqsort` by `HighsDomainChange::operator<`), so the 1.15.1
    /// consumers of `getImplications` keep their sorted input. The #3209 hash-tree merge is
    /// not ported.
    pub fn compute_implications(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: bool,
    ) -> bool {
        lend!(
            ms.mipdata_.implications,
            self,
            detach!(ms.mipdata_.domains[0], |gd| {
                gd.propagate(ms);
            })
        );
        if ms.mipdata_.domains[0].infeasible() || ms.mipdata_.domains[0].is_fixed(col) {
            return true;
        }

        // record redundant rows for lifting
        debug_assert!(
            ms.mipdata_.domains[0]
                .get_redundant_rows()
                .iter()
                .next()
                .is_none()
        );
        if self.store_lifting_opportunity.is_some() {
            ms.mipdata_.domains[0].set_record_redundant_rows(true);
        }

        let changedend: usize = ms.mipdata_.domains[0].get_changed_cols().len();

        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 30-33 at 8a5f9c1d8
        let dual_fix_probing_active: bool = ms.mipdata_.domains[0].get_dual_fix_probing_active();
        if dual_fix_probing_active {
            detach!(ms.mipdata_.domains[0], |gd| gd.dfp_begin_probing(ms));
        }

        let stackimplicstart: HighsInt =
            ms.mipdata_.domains[0].get_domain_change_stack().len() as HighsInt + 1;
        let mut numimplications_local: HighsInt = -stackimplicstart;
        if val {
            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.domains[0], |gd| {
                    gd.change_bound(
                        ms,
                        HighsDomainChange {
                            boundval: 1.0,
                            column: col,
                            boundtype: HighsBoundType::kLower,
                        },
                        HighsDomain_Reason::branching(),
                    );
                })
            );
        } else {
            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.domains[0], |gd| {
                    gd.change_bound(
                        ms,
                        HighsDomainChange {
                            boundval: 0.0,
                            column: col,
                            boundtype: HighsBoundType::kUpper,
                        },
                        HighsDomain_Reason::branching(),
                    );
                })
            );
        }

        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 54-57 at 8a5f9c1d8
        if dual_fix_probing_active && self.store_lifting_opportunity.is_some() {
            ms.mipdata_.domains[0]
                .get_dual_fix_probing_propagation()
                .store_lifting_opportunity = Some(std::sync::Arc::new(
                move |gd: &mut HighsDomain, ms: &mut HighsMipSolver| {
                    Self::store_lifting_opportunities_dfp(gd, ms, col, val)
                },
            ));
        }

        if self.is_infeasible_after_change(ms, col, val, changedend, dual_fix_probing_active) {
            return true;
        }

        lend!(
            ms.mipdata_.implications,
            self,
            detach!(ms.mipdata_.domains[0], |gd| {
                gd.propagate(ms);
            })
        );
        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 78-80 at 8a5f9c1d8
        if dual_fix_probing_active {
            ms.mipdata_.domains[0]
                .get_dual_fix_probing_propagation()
                .end_probing();
        }

        if self.is_infeasible_after_change(ms, col, val, changedend, dual_fix_probing_active) {
            return true;
        }

        let stackimplicend: HighsInt =
            ms.mipdata_.domains[0].get_domain_change_stack().len() as HighsInt;
        numimplications_local += stackimplicend;
        ms.mipdata_.pseudocosts[0].add_inference_observation(col, numimplications_local, val);

        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 89-209 at 8a5f9c1d8
        // (v1.16.0 lines 85-209). The new tail replaces the old one for every probe, not only
        // while dual-fix probing is active: implications and variable bounds are stored in
        // stack order and `implics` is sorted only at the end.
        if dual_fix_probing_active || ms.options_mip_.mip_dual_fix_probing {
            return self.compute_implications_dfp_tail(
                ms,
                col,
                val,
                (stackimplicstart, stackimplicend),
                numimplications_local,
                changedend,
                dual_fix_probing_active,
            );
        }

        let mut implics: Vec<HighsDomainChange> =
            Vec::with_capacity(numimplications_local.max(0) as usize);

        let num_entries: HighsInt = ms.mipdata_.cliquetable.get_num_entries();
        let max_entries: HighsInt = 100000 + ms.num_nonzero();

        for i in stackimplicstart..stackimplicend {
            let reason_type = ms.mipdata_.domains[0].get_domain_change_reason()[i as usize].r#type;
            let reason_index = ms.mipdata_.domains[0].get_domain_change_reason()[i as usize].index;
            if reason_type == HighsDomain_Reason::kCliqueTable
                && ((reason_index >> 1) == col || num_entries >= max_entries)
            {
                continue;
            }

            implics.push(ms.mipdata_.domains[0].get_domain_change_stack()[i as usize].clone());
        }

        // inform caller about lifting opportunities
        self.store_lifting_opportunities(ms, col, val);

        // backtrack
        self.do_backtrack(ms, changedend);

        // add the implications of binary variables to the clique table
        let binstart: usize = cpp_partition(&mut implics, |a| {
            !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, a.column))
        });

        pdqsort_by(&mut implics[..binstart], |a, b| a.op_lt(b));

        let mut clique: [HighsCliqueTable_CliqueVar; 2] = [
            HighsCliqueTable_CliqueVar::new(col, HighsInt::from(val)),
            HighsCliqueTable_CliqueVar::new(0, 0),
        ];

        for i in binstart..implics.len() {
            if implics[i].boundtype == HighsBoundType::kLower {
                clique[1] = HighsCliqueTable_CliqueVar::new(implics[i].column, 0);
            } else {
                clique[1] = HighsCliqueTable_CliqueVar::new(implics[i].column, 1);
            }

            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                    ms,
                    &mut clique,
                    false,
                    kHighsIInf
                ))
            );
            if ms.mipdata_.domains[0].infeasible() || ms.mipdata_.domains[0].is_fixed(col) {
                return true;
            }
        }

        // store variable bounds derived from implications
        for chg in implics[..binstart].iter() {
            let column = chg.column;
            let boundval = chg.boundval;
            if chg.boundtype == HighsBoundType::kLower {
                if val {
                    let lb = ms.mipdata_.domains[0].col_lower_[column as usize];
                    if lb != -kHighsInf {
                        self.add_vlb(ms, column, col, boundval - lb, lb);
                    }
                } else {
                    // in case the lower bound is infinite the varbound can
                    // still be tightened as soon as a finite upper bound is
                    // known because the offset is finite
                    let lb = ms.mipdata_.domains[0].col_lower_[column as usize];
                    self.add_vlb(ms, column, col, lb - boundval, boundval);
                }
            } else if val {
                let ub = ms.mipdata_.domains[0].col_upper_[column as usize];
                if ub != kHighsInf {
                    self.add_vub(ms, column, col, boundval - ub, ub);
                }
            } else {
                // in case the upper bound is infinite the varbound can
                // still be tightened as soon as a finite upper bound is
                // known because the offset is finite
                let ub = ms.mipdata_.domains[0].col_upper_[column as usize];
                self.add_vub(ms, column, col, ub - boundval, boundval);
            }
        }

        let loc: usize = (2 * col + HighsInt::from(val)) as usize;
        self.implications[loc].computed = true;
        implics.truncate(binstart);

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 143-154 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        if keep {
            self.has_probed[loc] = true;
            for chg in implics.iter() {
                let mut implication = HighsImplications_Implication::default();
                if chg.boundtype == HighsBoundType::kLower {
                    implication.lb = chg.boundval;
                } else {
                    implication.ub = chg.boundval;
                }
                self.add_implication(col, val, chg.column, implication);
            }
        }
        if !implics.is_empty() {
            let n = implics.len() as i64;
            self.implications[loc].implics = implics;
            if !keep {
                self.num_implications += n;
            }
        }

        false
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::computeImplications (lines 89-209)
    ///
    /// The latest tail of `computeImplications`, called when the global domain's dual-fix
    /// probing is active or `mip_dual_fix_probing` is on (`stack_range` =
    /// `(stackimplicstart, stackimplicend)`; `dual_fix_probing_active` is the C++
    /// `dualFixProbingActive`).
    ///
    /// PORT NOTE: `std::vector<HighsDomainChange>& implics = val ? implicationsUp :
    /// implicationsDown` is moved out of the member and moved back at every return, because
    /// the clique and variable-bound calls lend `self` back to `ms.mipdata_.implications`.
    #[allow(clippy::too_many_arguments)]
    fn compute_implications_dfp_tail(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: bool,
        stack_range: (HighsInt, HighsInt),
        numimplications_local: HighsInt,
        changedend: usize,
        dual_fix_probing_active: bool,
    ) -> bool {
        let (stackimplicstart, stackimplicend) = stack_range;
        // C++ `origStackCopy` is filled only `if (dualFixProbingActive)` (v1.16.0 lines 89-93)
        let orig_stack_copy: Vec<HighsDomainChange> = if dual_fix_probing_active {
            ms.mipdata_.domains[0].get_domain_change_stack()
                [stackimplicstart as usize..stackimplicend as usize]
                .to_vec()
        } else {
            Vec::new()
        };

        let mut implics: Vec<HighsDomainChange> = if val {
            std::mem::take(&mut self.implications_up)
        } else {
            std::mem::take(&mut self.implications_down)
        };
        implics.clear();
        implics.reserve(numimplications_local.max(0) as usize);

        let num_entries: HighsInt = ms.mipdata_.cliquetable.get_num_entries();
        let max_entries: HighsInt = 100000 + ms.num_nonzero();

        // C++ `dualFixProbingActive ? ...getZeroCostFixingPosition() : stackimplicend`
        let unsafe_stack_start: HighsInt = if dual_fix_probing_active {
            ms.mipdata_.domains[0]
                .get_dual_fix_probing_propagation()
                .get_zero_cost_fixing_position()
        } else {
            stackimplicend
        };
        let mut safe_implics_end: HighsInt = 0;

        for i in stackimplicstart..stackimplicend {
            let reason_type = ms.mipdata_.domains[0].get_domain_change_reason()[i as usize].r#type;
            let reason_index = ms.mipdata_.domains[0].get_domain_change_reason()[i as usize].index;
            if reason_type == HighsDomain_Reason::kCliqueTable
                && ((reason_index >> 1) == col || num_entries >= max_entries)
            {
                continue;
            }

            implics.push(ms.mipdata_.domains[0].get_domain_change_stack()[i as usize].clone());
            if i < unsafe_stack_start {
                safe_implics_end += 1;
            }
        }

        // inform caller about lifting opportunities
        self.store_lifting_opportunities(ms, col, val);

        // backtrack
        self.do_backtrack(ms, changedend);

        for implic in orig_stack_copy.iter() {
            if detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, implic.column)) {
                self.record_tentative_cliques(HighsInt::from(val), implic);
            }
        }

        if safe_implics_end < implics.len() as HighsInt {
            let safe = safe_implics_end as usize;
            let binstart: usize = safe
                + cpp_partition(&mut implics[safe..], |a| {
                    !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, a.column))
                });
            implics.truncate(binstart);
        }

        // add the implications of binary variables to the clique table
        let safe = safe_implics_end as usize;
        let binstart: usize = cpp_partition(&mut implics[..safe], |a| {
            !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, a.column))
        });

        let mut clique: [HighsCliqueTable_CliqueVar; 2] = [
            HighsCliqueTable_CliqueVar::new(col, HighsInt::from(val)),
            HighsCliqueTable_CliqueVar::new(0, 0),
        ];

        for i in binstart..safe {
            if implics[i].boundtype == HighsBoundType::kLower {
                clique[1] = HighsCliqueTable_CliqueVar::new(implics[i].column, 0);
            } else {
                clique[1] = HighsCliqueTable_CliqueVar::new(implics[i].column, 1);
            }

            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                    ms,
                    &mut clique,
                    false,
                    kHighsIInf
                ))
            );
            if ms.mipdata_.domains[0].infeasible() || ms.mipdata_.domains[0].is_fixed(col) {
                if val {
                    self.implications_up = implics;
                } else {
                    self.implications_down = implics;
                }
                return true;
            }
        }

        let num_erased_binaries: HighsInt = (safe - binstart) as HighsInt;
        implics.drain(binstart..safe);
        safe_implics_end -= num_erased_binaries;
        let safe = safe_implics_end as usize;

        // store variable bounds derived from implications
        for chg in implics[..safe].iter() {
            let column = chg.column;
            let boundval = chg.boundval;
            if chg.boundtype == HighsBoundType::kLower {
                if val {
                    let lb = ms.mipdata_.domains[0].col_lower_[column as usize];
                    if lb != -kHighsInf {
                        self.add_vlb(ms, column, col, boundval - lb, lb);
                    }
                } else {
                    // in case the lower bound is infinite the varbound can
                    // still be tightened as soon as a finite upper bound is
                    // known because the offset is finite
                    let lb = ms.mipdata_.domains[0].col_lower_[column as usize];
                    self.add_vlb(ms, column, col, lb - boundval, boundval);
                }
            } else if val {
                let ub = ms.mipdata_.domains[0].col_upper_[column as usize];
                if ub != kHighsInf {
                    self.add_vub(ms, column, col, boundval - ub, ub);
                }
            } else {
                // in case the upper bound is infinite the varbound can
                // still be tightened as soon as a finite upper bound is
                // known because the offset is finite
                let ub = ms.mipdata_.domains[0].col_upper_[column as usize];
                self.add_vub(ms, column, col, ub - boundval, boundval);
            }
        }

        // C++ `hasProbed[idx] = true` and `addImplication(idx, ...)` over the safe prefix
        // (lines 193-203), in the 1.15.1 storage (see compute_implications).
        let loc: usize = (2 * col + HighsInt::from(val)) as usize;
        self.implications[loc].computed = true;
        let mut stored: Vec<HighsDomainChange> = implics[..safe].to_vec();
        pdqsort_by(&mut stored, |a, b| a.op_lt(b));
        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 193-203 at 8a5f9c1d8
        // (the C++ loop walks the unsorted safe prefix `implics[..safe]`, not `stored`)
        let keep = ms.options_mip_.mip_keep_implications;
        if keep {
            self.has_probed[loc] = true;
            for chg in implics[..safe].iter() {
                let mut implication = HighsImplications_Implication::default();
                if chg.boundtype == HighsBoundType::kLower {
                    implication.lb = chg.boundval;
                } else {
                    implication.ub = chg.boundval;
                }
                self.add_implication(col, val, chg.column, implication);
            }
        }
        if !stored.is_empty() {
            let n = stored.len() as i64;
            self.implications[loc].implics = stored;
            if !keep {
                self.num_implications += n;
            }
        }

        pdqsort_by(&mut implics, |a, b| a.column < b.column);

        if val {
            self.implications_up = implics;
        } else {
            self.implications_down = implics;
        }

        false
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::computeImplications (lambda storeLiftingOpportunities via DualFixProbingPropagation::storeLiftingOpportunity, lines 42-57)
    ///
    /// PORT NOTE: the C++ closure `[&]() { storeLiftingOpportunities(col, val); }` is called
    /// by `propagateZeroCosts` inside `globaldomain.propagate()`. There the global domain is
    /// detached (`gd`) and this `HighsImplications` is lent back to
    /// `ms.mipdata_.implications` (the `propagate` calls in `compute_implications` run inside
    /// `lend!`), so the body reads the domain through `gd` and the callback buffer through
    /// `ms.mipdata_.implications`.
    fn store_lifting_opportunities_dfp(
        gd: &mut HighsDomain,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: bool,
    ) {
        if ms.mipdata_.implications.store_lifting_opportunity.is_some() {
            let rows: Vec<HighsInt> = gd
                .get_redundant_rows()
                .iter()
                .map(|elm| *elm.key())
                .collect();
            for row in rows {
                let value = gd.get_redundant_row_value(ms, row);
                let entry = (
                    row,
                    col,
                    if val { 1 } else { 0 },
                    (if val { -1.0 } else { 1.0 }) * value,
                );
                if let Some(buffer) = ms.mipdata_.implications.store_lifting_opportunity.as_mut() {
                    buffer.push(entry);
                }
            }
            gd.clear_redundant_rows();
            gd.set_record_redundant_rows(false);
        }
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::recordTentativeCliques (lines 272-278)
    fn record_tentative_cliques(&mut self, val: HighsInt, domchg: &HighsDomainChange) {
        let col = domchg.column;
        let fixing = &mut self.dual_fix_probing_bin_flags_[col as usize];
        if fixing.is_undecided() {
            self.dual_fix_probing_bin_inds_.push(col);
        }
        fixing.record(val == 1, domchg.boundtype);
    }

    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::clearTentativeCliques (lines 280-284)
    fn clear_tentative_cliques(&mut self) {
        for &col in self.dual_fix_probing_bin_inds_.iter() {
            self.dual_fix_probing_bin_flags_[col as usize].clear();
        }
        self.dual_fix_probing_bin_inds_.clear();
    }

    /// PORT NOTE (lambda `isInfeasible`, HighsImplications.cpp:53-60): returns true and
    /// restores the domain when the global domain became infeasible after a change.
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::computeImplications (lambda isInfeasible, lines 64-74)
    fn is_infeasible_after_change(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: bool,
        changedend: usize,
        dual_fix_probing_active: bool,
    ) -> bool {
        if !ms.mipdata_.domains[0].infeasible() {
            return false;
        }
        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 66-68 at 8a5f9c1d8
        if dual_fix_probing_active {
            ms.mipdata_.domains[0]
                .get_dual_fix_probing_propagation()
                .end_probing();
        }
        self.store_lifting_opportunities(ms, col, val);
        self.do_backtrack(ms, changedend);
        lend!(
            ms.mipdata_.implications,
            self,
            detach!(ms.mipdata_.cliquetable, |ct| {
                detach!(ms.mipdata_.domains[0], |gd| ct.vertex_infeasible(
                    ms,
                    gd,
                    col,
                    HighsInt::from(val)
                ));
            })
        );
        true
    }

    /// PORT NOTE (lambda `storeLiftingOpportunities`, HighsImplications.cpp:43-51).
    fn store_lifting_opportunities(&mut self, ms: &mut HighsMipSolver, col: HighsInt, val: bool) {
        if self.store_lifting_opportunity.is_some() {
            let rows: Vec<HighsInt> = ms.mipdata_.domains[0]
                .get_redundant_rows()
                .iter()
                .map(|elm| *elm.key())
                .collect();
            for row in rows {
                let value = detach!(ms.mipdata_.domains[0], |gd| gd
                    .get_redundant_row_value(ms, row));
                let entry = (
                    row,
                    col,
                    if val { 1 } else { 0 },
                    (if val { -1.0 } else { 1.0 }) * value,
                );
                if let Some(buffer) = self.store_lifting_opportunity.as_mut() {
                    buffer.push(entry);
                }
            }
            ms.mipdata_.domains[0].clear_redundant_rows();
            ms.mipdata_.domains[0].set_record_redundant_rows(false);
        }
    }

    /// PORT NOTE (lambda `doBacktrack`, HighsImplications.cpp:61-64).
    fn do_backtrack(&mut self, ms: &mut HighsMipSolver, changedend: usize) {
        detach!(ms.mipdata_.domains[0], |gd| {
            gd.backtrack(ms);
        });
        ms.mipdata_.domains[0].clear_changed_cols_2(changedend);
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::getBestVub (lines 157-228)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::getBestVub (lines 162-233)
    /// C++: `std::pair<HighsInt, VarBound> getBestVub(HighsInt col, const HighsSolution & lpSolution, double & bestUb, const HighsDomain & globaldom) const`
    pub fn get_best_vub(
        &self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        lp_solution: &HighsSolution,
        best_ub: &mut f64,
        globaldom: &HighsDomain,
    ) -> (HighsInt, HighsImplications_VarBound) {
        let feastol = ms.mipdata_.feastol;
        let inf_safe = ms.options_mip_.mip_aggregate_implications;
        let mut best_vub: (HighsInt, HighsImplications_VarBound) = (
            -1,
            HighsImplications_VarBound {
                coef: 0.0,
                constant: kHighsInf,
                origin: -1,
            },
        );
        let mut minbest_ub: f64 = *best_ub;
        let mut best_ub_dist: f64 = kHighsInf;
        let mut bestvubnodes: i64 = 0;

        let mut scale = globaldom.col_upper_[col as usize] - globaldom.col_lower_[col as usize];
        if scale == kHighsInf {
            scale = 1.0;
        } else {
            scale = 1.0 / scale;
        }

        // PORT NOTE (perf): the C++ `for_each` lambda, read-only, so
        // it runs on the tree directly (no snapshot); `return` is the lambda's `return`.
        let ms: &HighsMipSolver = ms;
        let col_value = &lp_solution.col_value;
        self.vubs[col as usize].for_each(&mut |&vub_col, vub| {
            let vub = *vub;
            if vub.coef == kHighsInf {
                return;
            }
            if globaldom.is_fixed(vub_col) {
                return;
            }
            debug_assert!(globaldom.is_binary(ms, vub_col));
            let vub_val_col = col_value[vub_col as usize];
            let vubval = vub_val_col * vub.coef + vub.constant;
            let mut ub_dist = cpp_max(0.0, vubval - col_value[col as usize]);

            let y_dist = feastol
                + (if vub.coef > 0.0 {
                    1.0 - vub_val_col
                } else {
                    vub_val_col
                });

            // skip variable bound if the distance towards the variable bound constraint
            // is larger than the distance to the point where the binary column is
            // relaxed to it's weakest bound, i.e. 1 if its coefficient is positive.
            // The variable bound constraint has the form x <= ay + b with y binary and
            // hence the norm is sqrt(1 + a^2) and the distance of the variable bound
            // constraint is ay + b - x evaluated at the solution values of x and y
            // divided by the norm.
            let norm2 = 1.0 + vub.coef * vub.coef;
            if K_SKIP_BAD_VBDS && ub_dist * ub_dist > y_dist * y_dist * norm2 {
                return;
            }

            debug_assert!(vub_col >= 0 && vub_col < ms.num_col());
            ub_dist *= scale;
            if ub_dist <= best_ub_dist + feastol {
                let minvubval = var_bound_min_value_x(&vub, inf_safe);
                let vubnodes: i64 = if vub.coef > 0.0 {
                    ms.mipdata_.nodequeue.num_nodes_down(vub_col)
                } else {
                    ms.mipdata_.nodequeue.num_nodes_up(vub_col)
                };

                if is_vub_better(
                    ub_dist,
                    vubnodes,
                    minvubval,
                    vub_col,
                    &vub,
                    &best_vub,
                    best_ub_dist,
                    bestvubnodes,
                    minbest_ub,
                    lp_solution,
                    feastol,
                ) {
                    *best_ub = vubval;
                    minbest_ub = minvubval;
                    best_vub = (vub_col, vub);
                    bestvubnodes = vubnodes;
                    best_ub_dist = ub_dist;
                }
            }
        });

        best_vub
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::getBestVlb (lines 230-296)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::getBestVlb (lines 235-301)
    /// C++: `std::pair<HighsInt, VarBound> getBestVlb(HighsInt col, const HighsSolution & lpSolution, double & bestLb, const HighsDomain & globaldom) const`
    pub fn get_best_vlb(
        &self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        lp_solution: &HighsSolution,
        best_lb: &mut f64,
        globaldom: &HighsDomain,
    ) -> (HighsInt, HighsImplications_VarBound) {
        let feastol = ms.mipdata_.feastol;
        let inf_safe = ms.options_mip_.mip_aggregate_implications;
        let mut best_vlb: (HighsInt, HighsImplications_VarBound) = (
            -1,
            HighsImplications_VarBound {
                coef: 0.0,
                constant: -kHighsInf,
                origin: -1,
            },
        );
        let mut maxbest_lb: f64 = *best_lb;
        let mut best_lb_dist: f64 = kHighsInf;
        let mut bestvlbnodes: i64 = 0;

        let mut scale = globaldom.col_upper_[col as usize] - globaldom.col_lower_[col as usize];
        if scale == kHighsInf {
            scale = 1.0;
        } else {
            scale = 1.0 / scale;
        }

        // PORT NOTE (perf): read-only `for_each`, as in `get_best_vub`.
        let ms: &HighsMipSolver = ms;
        let col_value = &lp_solution.col_value;
        self.vlbs[col as usize].for_each(&mut |&vlb_col, vlb| {
            let vlb = *vlb;
            if vlb.coef == -kHighsInf {
                return;
            }
            if globaldom.is_fixed(vlb_col) {
                return;
            }
            debug_assert!(globaldom.is_binary(ms, vlb_col));
            debug_assert!(vlb_col >= 0 && vlb_col < ms.num_col());
            let vlb_val_col = col_value[vlb_col as usize];
            let vlbval = vlb_val_col * vlb.coef + vlb.constant;
            let mut lb_dist = cpp_max(0.0, col_value[col as usize] - vlbval);

            let y_dist = feastol
                + (if vlb.coef > 0.0 {
                    vlb_val_col
                } else {
                    1.0 - vlb_val_col
                });

            let norm2 = 1.0 + vlb.coef * vlb.coef;
            if K_SKIP_BAD_VBDS && lb_dist * lb_dist > y_dist * y_dist * norm2 {
                return;
            }

            // scale the distance as if the bounded column was scaled to have ub-lb=1
            lb_dist *= scale;
            if lb_dist <= best_lb_dist + feastol {
                let maxvlbval = var_bound_max_value_x(&vlb, inf_safe);
                let vlbnodes: i64 = if vlb.coef > 0.0 {
                    ms.mipdata_.nodequeue.num_nodes_up(vlb_col)
                } else {
                    ms.mipdata_.nodequeue.num_nodes_down(vlb_col)
                };

                if is_vlb_better(
                    lb_dist,
                    vlbnodes,
                    maxvlbval,
                    vlb_col,
                    &vlb,
                    &best_vlb,
                    best_lb_dist,
                    bestvlbnodes,
                    maxbest_lb,
                    lp_solution,
                    feastol,
                ) {
                    *best_lb = vlbval;
                    maxbest_lb = maxvlbval;
                    best_vlb = (vlb_col, vlb);
                    bestvlbnodes = vlbnodes;
                    best_lb_dist = lb_dist;
                }
            }
        });

        best_vlb
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::runProbing (lines 298-389)
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::runProbing (lines 357-495)
    /// C++: `bool runProbing(HighsInt col, HighsInt & numReductions)`
    ///
    /// PORT NOTE: `getImplications(col, 0)` and `getImplications(col, 1)` each borrow `self`
    /// mutably in Rust, so both lists are copied (`to_vec`) before the merge. Both are cached
    /// at this point, so the copy does not change any result.
    ///
    /// PORT NOTE (HiGHS #3244): with the global domain's dual-fix probing active, the merge reads
    /// the latest `implicationsDown` / `implicationsUp` (copied likewise; the merge does not
    /// write them) and returns once the domain is infeasible after a bound change. The latest
    /// `probedBefore` is the 1.15.1 `implicationsCached` in the kept storage, and the closing
    /// `hasProbed[ImplIdx{col, 0/1}] = true` (lines 490-491) is the `computed` flag that
    /// `compute_implications` has already set for both values.
    pub fn run_probing(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        num_reductions: &mut HighsInt,
    ) -> bool {
        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 305-307 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        let col_is_binary = detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, col));
        if col_is_binary
            && !(if keep {
                self.probed_before(col, true)
            } else {
                self.implications_cached(col, true)
            })
            && !(if keep {
                self.probed_before(col, false)
            } else {
                self.implications_cached(col, false)
            })
            && ms.mipdata_.cliquetable.get_substitution(col).is_none()
        {
            // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 362-367 at 8a5f9c1d8
            let dual_fix_probing_active: bool =
                ms.mipdata_.domains[0].get_dual_fix_probing_active();
            if dual_fix_probing_active {
                self.clear_tentative_cliques();
                ms.mipdata_.domains[0]
                    .get_dual_fix_probing_propagation()
                    .set_zero_cost_fixing_position(kHighsIInf);
            }

            let mut infeasible = self.compute_implications(ms, col, true);
            if ms.mipdata_.domains[0].infeasible() {
                return true;
            }
            if infeasible {
                return true;
            }
            if ms.mipdata_.cliquetable.get_substitution(col).is_some() {
                return true;
            }

            infeasible = self.compute_implications(ms, col, false);
            if ms.mipdata_.domains[0].infeasible() {
                return true;
            }
            if infeasible {
                return true;
            }
            if ms.mipdata_.cliquetable.get_substitution(col).is_some() {
                return true;
            }

            // HiGHS #3244 (mip_dual_fix_probing): v1.16.0 runProbing (lines 419-485) reads the
            // `implicationsDown`/`implicationsUp` members for every probe; the f7b87ae01 tree walk
            // below is only the pre-#3244 state.
            let x1_dfp = ms.options_mip_.mip_dual_fix_probing;
            // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 320-373 at f7b87ae01
            if keep && !dual_fix_probing_active && !x1_dfp {
                let mut entries: Vec<(HighsInt, HighsImplications_Implication)> = Vec::new();
                self.implication_tree[(2 * col) as usize]
                    .for_each(&mut |k, v| entries.push((*k, *v)));
                for (implcol, down_impl) in entries.into_iter() {
                    if self.colsubstituted[implcol as usize] != 0
                        || ms.mipdata_.domains[0].is_fixed(implcol)
                    {
                        continue;
                    }
                    let up = self.implication_tree[(2 * col + 1) as usize]
                        .find(&implcol)
                        .copied();
                    if let Some(up_impl) = up {
                        let feastol = ms.mipdata_.feastol;
                        let lb_down = cpp_max(
                            ms.mipdata_.domains[0].col_lower_[implcol as usize],
                            down_impl.lb,
                        );
                        let ub_down = cpp_min(
                            ms.mipdata_.domains[0].col_upper_[implcol as usize],
                            down_impl.ub,
                        );
                        let lb_up = cpp_max(
                            ms.mipdata_.domains[0].col_lower_[implcol as usize],
                            up_impl.lb,
                        );
                        let ub_up = cpp_min(
                            ms.mipdata_.domains[0].col_upper_[implcol as usize],
                            up_impl.ub,
                        );

                        if lb_down == ub_down && lb_up == ub_up && (lb_down - lb_up).abs() > feastol
                        {
                            let substitution = HighsSubstitution {
                                substcol: implcol,
                                staycol: col,
                                offset: lb_down,
                                scale: lb_up - lb_down,
                            };
                            self.substitutions.push(substitution);
                            self.colsubstituted[implcol as usize] = 1;
                            *num_reductions += 1;
                        } else if !ms.parallel_lock_active() {
                            let lb = cpp_min(lb_down, lb_up);
                            let ub = cpp_max(ub_down, ub_up);

                            if lb > ms.mipdata_.domains[0].col_lower_[implcol as usize] {
                                lend!(
                                    ms.mipdata_.implications,
                                    self,
                                    detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                                        ms,
                                        HighsDomainChange {
                                            boundval: lb,
                                            column: implcol,
                                            boundtype: HighsBoundType::kLower,
                                        },
                                        HighsDomain_Reason::unspecified(),
                                    ))
                                );
                                *num_reductions += 1;
                            }

                            if ub < ms.mipdata_.domains[0].col_upper_[implcol as usize] {
                                lend!(
                                    ms.mipdata_.implications,
                                    self,
                                    detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                                        ms,
                                        HighsDomainChange {
                                            boundval: ub,
                                            column: implcol,
                                            boundtype: HighsBoundType::kUpper,
                                        },
                                        HighsDomain_Reason::unspecified(),
                                    ))
                                );
                                *num_reductions += 1;
                            }
                        }
                    }
                    if ms.mipdata_.domains[0].infeasible() {
                        break;
                    }
                }
                self.has_probed[(2 * col) as usize] = true;
                self.has_probed[(2 * col + 1) as usize] = true;
                return true;
            }

            let implicsdown: Vec<HighsDomainChange>;
            let implicsup: Vec<HighsDomainChange>;
            // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 381-420 at 8a5f9c1d8
            if dual_fix_probing_active {
                if !self.dual_fix_probing_bin_inds_.is_empty() {
                    let mut clique: [HighsCliqueTable_CliqueVar; 2] = [
                        HighsCliqueTable_CliqueVar::new(0, 0),
                        HighsCliqueTable_CliqueVar::new(0, 0),
                    ];
                    // C++ range-for over dualFixProbingBinInds_; the body never pushes to it.
                    for ki in 0..self.dual_fix_probing_bin_inds_.len() {
                        let k = self.dual_fix_probing_bin_inds_[ki];
                        if !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, k))
                            || self.colsubstituted[k as usize] != 0
                        {
                            continue;
                        }
                        if ms.mipdata_.domains[0].infeasible() {
                            return true;
                        }
                        let f = self.dual_fix_probing_bin_flags_[k as usize];
                        if f.is_undecided() {
                            continue;
                        }
                        if f.down_probe == HighsImplications_TentativeFixingDirection::FixLower
                            && f.up_probe == HighsImplications_TentativeFixingDirection::FixLower
                        {
                            let v = ms.mipdata_.domains[0].col_lower_[k as usize];
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.domains[0], |gd| gd.fix_col(
                                    ms,
                                    k,
                                    v,
                                    HighsDomain_Reason::unspecified()
                                ))
                            );
                        } else if f.down_probe
                            == HighsImplications_TentativeFixingDirection::FixUpper
                            && f.up_probe == HighsImplications_TentativeFixingDirection::FixUpper
                        {
                            let v = ms.mipdata_.domains[0].col_upper_[k as usize];
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.domains[0], |gd| gd.fix_col(
                                    ms,
                                    k,
                                    v,
                                    HighsDomain_Reason::unspecified()
                                ))
                            );
                        } else if f.down_probe
                            == HighsImplications_TentativeFixingDirection::FixUpper
                            && f.up_probe == HighsImplications_TentativeFixingDirection::FixLower
                            && !ms.mipdata_.cliquetable.is_full()
                        {
                            clique[0] = HighsCliqueTable_CliqueVar::new(col, 1);
                            clique[1] = HighsCliqueTable_CliqueVar::new(k, 1);
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                                    ms,
                                    &mut clique,
                                    false,
                                    kHighsIInf
                                ))
                            );
                            clique[0] = HighsCliqueTable_CliqueVar::new(col, 0);
                            clique[1] = HighsCliqueTable_CliqueVar::new(k, 0);
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                                    ms,
                                    &mut clique,
                                    false,
                                    kHighsIInf
                                ))
                            );
                        } else if f.down_probe
                            == HighsImplications_TentativeFixingDirection::FixLower
                            && f.up_probe == HighsImplications_TentativeFixingDirection::FixUpper
                            && !ms.mipdata_.cliquetable.is_full()
                        {
                            clique[0] = HighsCliqueTable_CliqueVar::new(col, 1);
                            clique[1] = HighsCliqueTable_CliqueVar::new(k, 0);
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                                    ms,
                                    &mut clique,
                                    false,
                                    kHighsIInf
                                ))
                            );
                            clique[0] = HighsCliqueTable_CliqueVar::new(col, 0);
                            clique[1] = HighsCliqueTable_CliqueVar::new(k, 1);
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.cliquetable, |ct| ct.add_clique(
                                    ms,
                                    &mut clique,
                                    false,
                                    kHighsIInf
                                ))
                            );
                        }
                        if ms.mipdata_.domains[0].infeasible() {
                            return true;
                        }
                    }

                    self.clear_tentative_cliques();
                }

                implicsdown = self.implications_down.to_vec();
                implicsup = self.implications_up.to_vec();
            } else if x1_dfp {
                implicsdown = self.implications_down.to_vec();
                implicsup = self.implications_up.to_vec();
            } else {
                // analyze implications
                implicsdown = self
                    .get_implications(ms, col, false, &mut infeasible)
                    .to_vec();
                implicsup = self
                    .get_implications(ms, col, true, &mut infeasible)
                    .to_vec();
            }
            let nimplicsdown = implicsdown.len();
            let nimplicsup = implicsup.len();
            let mut u: usize = 0;
            let mut d: usize = 0;

            while u < nimplicsup && d < nimplicsdown {
                if implicsup[u].column < implicsdown[d].column {
                    u += 1;
                } else if implicsdown[d].column < implicsup[u].column {
                    d += 1;
                } else {
                    let implcol = implicsup[u].column;
                    let mut lb_down = ms.mipdata_.domains[0].col_lower_[implcol as usize];
                    let mut ub_down = ms.mipdata_.domains[0].col_upper_[implcol as usize];
                    let mut lb_up = lb_down;
                    let mut ub_up = ub_down;

                    loop {
                        if implicsdown[d].boundtype == HighsBoundType::kLower {
                            lb_down = cpp_max(lb_down, implicsdown[d].boundval);
                        } else {
                            ub_down = cpp_min(ub_down, implicsdown[d].boundval);
                        }
                        d += 1;
                        if !(d < nimplicsdown && implicsdown[d].column == implcol) {
                            break;
                        }
                    }

                    loop {
                        if implicsup[u].boundtype == HighsBoundType::kLower {
                            lb_up = cpp_max(lb_up, implicsup[u].boundval);
                        } else {
                            ub_up = cpp_min(ub_up, implicsup[u].boundval);
                        }
                        u += 1;
                        if !(u < nimplicsup && implicsup[u].column == implcol) {
                            break;
                        }
                    }

                    if self.colsubstituted[implcol as usize] != 0
                        || ms.mipdata_.domains[0].is_fixed(implcol)
                    {
                        continue;
                    }

                    let feastol = ms.mipdata_.feastol;
                    if lb_down == ub_down && lb_up == ub_up && (lb_down - lb_up).abs() > feastol {
                        let substitution = HighsSubstitution {
                            substcol: implcol,
                            staycol: col,
                            offset: lb_down,
                            scale: lb_up - lb_down,
                        };
                        self.substitutions.push(substitution);
                        self.colsubstituted[implcol as usize] = 1;
                        *num_reductions += 1;
                    } else if !ms.parallel_lock_active() {
                        let lb = cpp_min(lb_down, lb_up);
                        let ub = cpp_max(ub_down, ub_up);

                        if lb > ms.mipdata_.domains[0].col_lower_[implcol as usize] {
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                                    ms,
                                    HighsDomainChange {
                                        boundval: lb,
                                        column: implcol,
                                        boundtype: HighsBoundType::kLower,
                                    },
                                    HighsDomain_Reason::unspecified(),
                                ))
                            );
                            *num_reductions += 1;
                            // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp line 476 at 8a5f9c1d8
                            if (dual_fix_probing_active || x1_dfp)
                                && ms.mipdata_.domains[0].infeasible()
                            {
                                return true;
                            }
                        }

                        if ub < ms.mipdata_.domains[0].col_upper_[implcol as usize] {
                            lend!(
                                ms.mipdata_.implications,
                                self,
                                detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                                    ms,
                                    HighsDomainChange {
                                        boundval: ub,
                                        column: implcol,
                                        boundtype: HighsBoundType::kUpper,
                                    },
                                    HighsDomain_Reason::unspecified(),
                                ))
                            );
                            *num_reductions += 1;
                            // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp line 483 at 8a5f9c1d8
                            if (dual_fix_probing_active || x1_dfp)
                                && ms.mipdata_.domains[0].infeasible()
                            {
                                return true;
                            }
                        }
                    }
                }
            }

            return true;
        }

        false
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::strengthenVarBound (lines 391-409)
    /// C++: `void strengthenVarBound(VarBound & vbnd, HighsInt multiplier) const`
    pub fn strengthen_var_bound(
        &self,
        vbnd: &mut HighsImplications_VarBound,
        multiplier: HighsInt,
    ) {
        // try to strengthen variable bound constraint x + a * y <= b by computing
        // MIR cut. since x is a general-integer variable (with integral bounds) and
        // its coefficient is 1.0, it is not shifted (or complemented). similarly,
        // since y is a binary variable (i.e., its lower bound is 0.0), it does not
        // need to be shifted, and we also do not try to complement it.
        if vbnd.coef.abs() == kHighsInf || vbnd.constant.abs() == kHighsInf {
            return;
        }
        const F0MIN: f64 = 0.005;
        const F0MAX: f64 = 0.995;
        let m = multiplier as f64;
        let downrhs = (m * vbnd.constant).floor();
        let f0 = m * vbnd.constant - downrhs;
        if f0 < F0MIN || f0 > F0MAX {
            return;
        }
        let downaj = (((-multiplier) as f64) * vbnd.coef + kHighsTiny).floor();
        let fj = ((-multiplier) as f64) * vbnd.coef - downaj;
        vbnd.constant = m * downrhs;
        vbnd.coef = ((-multiplier) as f64) * (downaj + cpp_max(fj - f0, 0.0) / (1.0 - f0));
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::addVUB (lines 411-416)
    /// C++: `void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant)`
    pub fn add_vub(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vubcol: HighsInt,
        vubcoef: f64,
        vubconstant: f64,
    ) {
        let colupperbound = ms.mipdata_.domains[0].col_upper_[col as usize];
        let colisintegral = ms.is_col_integral(col);
        self.add_vub_2(
            ms,
            col,
            vubcol,
            vubcoef,
            vubconstant,
            colupperbound,
            colisintegral,
        );
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::addVUB (lines 418-450)
    /// C++: `void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant, double colupperbound, bool colisintegral)`
    pub fn add_vub_2(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vubcol: HighsInt,
        vubcoef: f64,
        vubconstant: f64,
        colupperbound: f64,
        colisintegral: bool,
    ) {
        // HiGHS #3264 (mip_aggregate_implications): HighsImplications.cpp lines 418-430 at
        // f2ae2b320 fold this overload into the `origin` one (default -1), with the inf-safe
        // `VarBound` values.
        if ms.options_mip_.mip_aggregate_implications {
            self.add_vub_2_o(
                ms,
                col,
                vubcol,
                vubcoef,
                vubconstant,
                colupperbound,
                colisintegral,
                -1,
            );
            return;
        }
        // assume that VUBs do not have infinite coefficients and infinite constant
        // terms since such VUBs effectively evaluate to NaN.
        debug_assert!(vubcoef.abs() != kHighsInf || vubconstant.abs() != kHighsInf);
        if self.too_many_var_bounds() {
            return;
        }

        let mut vub = HighsImplications_VarBound {
            coef: vubcoef,
            constant: vubconstant,
            origin: -1,
        };

        if colisintegral {
            // try to strengthen VUB
            self.strengthen_var_bound(&mut vub, 1);
            if vub.coef == 0.0 {
                return;
            }
        }

        ms.mipdata_
            .debug_solution
            .check_vub(col, vubcol, vubcoef, vubconstant);

        let feastol = ms.mipdata_.feastol;
        let min_bound = var_bound_min_value(&vub);
        if min_bound >= colupperbound - feastol {
            return;
        }

        let (currentvub, inserted) =
            self.vubs[col as usize].insert_or_get(HighsHashTableEntry::new(vubcol, vub));

        if !inserted {
            let current_min_bound = var_bound_min_value(currentvub);
            if min_bound < current_min_bound - feastol {
                currentvub.coef = vub.coef;
                currentvub.constant = vub.constant;
            }
        } else {
            self.num_var_bounds += 1;
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::addVUB (lines 432-467)
    /// C++: `void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant, double colupperbound, bool colisintegral, HighsInt origin)`
    ///
    /// PORT NOTE: the 1.15.1 `add_vub` / `add_vub_2` stay as they are (rule F): they use the old
    /// min/max value, which is NaN for an infinite coefficient (reached from
    /// `compute_implications`). Only this origin overload uses the latest inf-safe one.
    pub fn add_vub_2_o(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vubcol: HighsInt,
        vubcoef: f64,
        vubconstant: f64,
        colupperbound: f64,
        colisintegral: bool,
        origin: HighsInt,
    ) {
        // assume that VUBs do not have infinite coefficients and infinite constant
        // terms since such VUBs effectively evaluate to NaN.
        debug_assert!(vubcoef.abs() != kHighsInf || vubconstant.abs() != kHighsInf);
        if self.too_many_var_bounds() {
            return;
        }

        let mut vub = HighsImplications_VarBound {
            coef: vubcoef,
            constant: vubconstant,
            origin,
        };

        if colisintegral {
            // try to strengthen VUB
            self.strengthen_var_bound(&mut vub, 1);
            if vub.coef == 0.0 {
                return;
            }
        }

        ms.mipdata_
            .debug_solution
            .check_vub(col, vubcol, vubcoef, vubconstant);

        let feastol = ms.mipdata_.feastol;
        let min_bound = var_bound_min_value_inf_safe(&vub);
        if min_bound >= colupperbound - feastol {
            return;
        }

        let (currentvub, inserted) =
            self.vubs[col as usize].insert_or_get(HighsHashTableEntry::new(vubcol, vub));

        if !inserted {
            let current_min_bound = var_bound_min_value_inf_safe(currentvub);
            if min_bound < current_min_bound - feastol {
                currentvub.coef = vub.coef;
                currentvub.constant = vub.constant;
                currentvub.origin = vub.origin;
            }
        } else {
            self.num_var_bounds += 1;
        }

        if origin >= 0 {
            self.row_to_var_bounds[origin as usize].insert(HighsHashTableEntry::new(col, vubcol));
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::addVUB (lines 425-430)
    /// C++: `void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant, HighsInt origin)`
    pub fn add_vub_o(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vubcol: HighsInt,
        vubcoef: f64,
        vubconstant: f64,
        origin: HighsInt,
    ) {
        let colupperbound = ms.mipdata_.domains[0].col_upper_[col as usize];
        let colisintegral = ms.is_col_integral(col);
        self.add_vub_2_o(
            ms,
            col,
            vubcol,
            vubcoef,
            vubconstant,
            colupperbound,
            colisintegral,
            origin,
        );
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::addVLB (lines 452-457)
    /// C++: `void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant)`
    pub fn add_vlb(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vlbcol: HighsInt,
        vlbcoef: f64,
        vlbconstant: f64,
    ) {
        let colllowerbound = ms.mipdata_.domains[0].col_lower_[col as usize];
        let colisintegral = ms.is_col_integral(col);
        self.add_vlb_2(
            ms,
            col,
            vlbcol,
            vlbcoef,
            vlbconstant,
            colllowerbound,
            colisintegral,
        );
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::addVLB (lines 459-492)
    /// C++: `void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant, double colllowerbound, bool colisintegral)`
    pub fn add_vlb_2(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vlbcol: HighsInt,
        vlbcoef: f64,
        vlbconstant: f64,
        colllowerbound: f64,
        colisintegral: bool,
    ) {
        // HiGHS #3264 (mip_aggregate_implications): HighsImplications.cpp lines 418-430 at
        // f2ae2b320 fold this overload into the `origin` one (default -1), with the inf-safe
        // `VarBound` values.
        if ms.options_mip_.mip_aggregate_implications {
            self.add_vlb_2_o(
                ms,
                col,
                vlbcol,
                vlbcoef,
                vlbconstant,
                colllowerbound,
                colisintegral,
                -1,
            );
            return;
        }
        // assume that VLBs do not have infinite coefficients and infinite constant
        // terms since such VLBs effectively evaluate to NaN.
        debug_assert!(vlbcoef.abs() != kHighsInf || vlbconstant.abs() != kHighsInf);
        if self.too_many_var_bounds() {
            return;
        }

        let mut vlb = HighsImplications_VarBound {
            coef: vlbcoef,
            constant: vlbconstant,
            origin: -1,
        };

        if colisintegral {
            // try to strengthen VLB
            self.strengthen_var_bound(&mut vlb, -1);
            if vlb.coef == 0.0 {
                return;
            }
        }

        ms.mipdata_
            .debug_solution
            .check_vlb(col, vlbcol, vlbcoef, vlbconstant);

        let feastol = ms.mipdata_.feastol;
        let max_bound = var_bound_max_value(&vlb);
        if max_bound <= colllowerbound + feastol {
            return;
        }

        let (currentvlb, inserted) =
            self.vlbs[col as usize].insert_or_get(HighsHashTableEntry::new(vlbcol, vlb));

        if !inserted {
            let current_max_bound = var_bound_max_value(currentvlb);
            if max_bound > current_max_bound + feastol {
                currentvlb.coef = vlb.coef;
                currentvlb.constant = vlb.constant;
            }
        } else {
            self.num_var_bounds += 1;
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::addVLB (lines 476-512)
    /// C++: `void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant, double colllowerbound, bool colisintegral, HighsInt origin)`
    ///
    /// PORT NOTE: the 1.15.1 `add_vlb` / `add_vlb_2` stay as they are (rule F): they use the old
    /// min/max value, which is NaN for an infinite coefficient (reached from
    /// `compute_implications`). Only this origin overload uses the latest inf-safe one.
    pub fn add_vlb_2_o(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vlbcol: HighsInt,
        vlbcoef: f64,
        vlbconstant: f64,
        colllowerbound: f64,
        colisintegral: bool,
        origin: HighsInt,
    ) {
        // assume that VLBs do not have infinite coefficients and infinite constant
        // terms since such VLBs effectively evaluate to NaN.
        debug_assert!(vlbcoef.abs() != kHighsInf || vlbconstant.abs() != kHighsInf);
        if self.too_many_var_bounds() {
            return;
        }

        let mut vlb = HighsImplications_VarBound {
            coef: vlbcoef,
            constant: vlbconstant,
            origin,
        };

        if colisintegral {
            // try to strengthen VLB
            self.strengthen_var_bound(&mut vlb, -1);
            if vlb.coef == 0.0 {
                return;
            }
        }

        ms.mipdata_
            .debug_solution
            .check_vlb(col, vlbcol, vlbcoef, vlbconstant);

        let feastol = ms.mipdata_.feastol;
        let max_bound = var_bound_max_value_inf_safe(&vlb);
        if max_bound <= colllowerbound + feastol {
            return;
        }

        let (currentvlb, inserted) =
            self.vlbs[col as usize].insert_or_get(HighsHashTableEntry::new(vlbcol, vlb));

        if !inserted {
            let current_max_bound = var_bound_max_value_inf_safe(currentvlb);
            if max_bound > current_max_bound + feastol {
                currentvlb.coef = vlb.coef;
                currentvlb.constant = vlb.constant;
                currentvlb.origin = vlb.origin;
            }
        } else {
            self.num_var_bounds += 1;
        }

        if origin >= 0 {
            self.row_to_var_bounds[origin as usize].insert(HighsHashTableEntry::new(col, vlbcol));
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::addVLB (lines 469-474)
    /// C++: `void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant, HighsInt origin)`
    pub fn add_vlb_o(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vlbcol: HighsInt,
        vlbcoef: f64,
        vlbconstant: f64,
        origin: HighsInt,
    ) {
        let colllowerbound = ms.mipdata_.domains[0].col_lower_[col as usize];
        let colisintegral = ms.is_col_integral(col);
        self.add_vlb_2_o(
            ms,
            col,
            vlbcol,
            vlbcoef,
            vlbconstant,
            colllowerbound,
            colisintegral,
            origin,
        );
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::rebuild (lines 494-559)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::rebuild (lines 535-628)
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.cpp::HighsImplications::rebuild (lines 656-752)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::rebuild (lines 498-586)
    /// C++: `void rebuild(HighsInt ncols, const std::vector<HighsInt> & orig2reducedcol, const std::vector<HighsInt> & orig2reducedrow)`
    pub fn rebuild(
        &mut self,
        ms: &mut HighsMipSolver,
        ncols: HighsInt,
        orig2reducedcol: &[HighsInt],
        orig2reducedrow: &[HighsInt],
    ) {
        let oldvlbs = std::mem::take(&mut self.vlbs);
        let oldvubs = std::mem::take(&mut self.vubs);

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 507-508 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        let old_implication_tree = if keep {
            std::mem::take(&mut self.implication_tree)
        } else {
            Vec::new()
        };

        self.colsubstituted.clear();
        self.colsubstituted.shrink_to_fit();
        self.implications.clear();
        self.implications.shrink_to_fit();

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 514-517 at f7b87ae01
        if keep {
            self.has_probed.clear();
            self.has_probed.shrink_to_fit();
            self.reverse_implications.clear();
            self.reverse_implications.shrink_to_fit();
        }

        let ncols_us = ncols as usize;
        self.implications
            .resize_with(2 * ncols_us, HighsImplications_Implics::default);

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 519-521 at f7b87ae01
        if keep {
            self.implication_tree
                .resize_with(2 * ncols_us, HighsHashTree::default);
            self.has_probed.resize(2 * ncols_us, false);
            self.reverse_implications
                .resize_with(ncols_us, HighsHashTree::default);
        }

        self.colsubstituted.resize(ncols_us, 0);
        self.substitutions.clear();
        self.vubs.clear();
        self.vubs.shrink_to_fit();
        self.vubs.resize_with(ncols_us, HighsHashTree::default);
        self.vlbs.clear();
        self.vlbs.shrink_to_fit();
        self.vlbs.resize_with(ncols_us, HighsHashTree::default);
        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.cpp lines 688-690 at 8a5f9c1d8
        if ms.options_mip_.mip_dual_fix_probing {
            self.dual_fix_probing_bin_inds_.clear();
            self.dual_fix_probing_bin_inds_.reserve(ncols_us);
            self.dual_fix_probing_bin_flags_.clear();
            self.dual_fix_probing_bin_flags_
                .resize(ncols_us, HighsImplications_TentativeFixing::default());
        }
        self.row_to_var_bounds.clear();
        self.row_to_var_bounds.shrink_to_fit();
        self.resize_rows(ms.num_row());
        self.num_implications = 0;
        self.num_var_bounds = 0;

        self.next_cleanup_call = ms.num_nonzero();

        for (i, (oldvub_i, oldvlb_i)) in oldvubs.iter().zip(oldvlbs.iter()).enumerate() {
            let newi = orig2reducedcol[i];

            if newi == -1
                || !ms
                    .mipdata_
                    .post_solve_stack
                    .is_col_linearly_transformable(newi)
            {
                continue;
            }

            oldvub_i.for_each(&mut |&vub_col, vub| {
                let vub = *vub;
                let new_vub_col = orig2reducedcol[vub_col as usize];
                if new_vub_col == -1 {
                    return;
                }

                if !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, new_vub_col))
                    || !ms
                        .mipdata_
                        .post_solve_stack
                        .is_col_linearly_transformable(new_vub_col)
                {
                    return;
                }

                if vub.origin >= 0 {
                    self.add_vub_o(
                        ms,
                        newi,
                        new_vub_col,
                        vub.coef,
                        vub.constant,
                        orig2reducedrow[vub.origin as usize],
                    );
                } else {
                    self.add_vub(ms, newi, new_vub_col, vub.coef, vub.constant);
                }
            });

            oldvlb_i.for_each(&mut |&vlb_col, vlb| {
                let vlb = *vlb;
                let new_vlb_col = orig2reducedcol[vlb_col as usize];
                if new_vlb_col == -1 {
                    return;
                }

                if !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, new_vlb_col))
                    || !ms
                        .mipdata_
                        .post_solve_stack
                        .is_col_linearly_transformable(new_vlb_col)
                {
                    return;
                }

                if vlb.origin >= 0 {
                    self.add_vlb_o(
                        ms,
                        newi,
                        new_vlb_col,
                        vlb.coef,
                        vlb.constant,
                        orig2reducedrow[vlb.origin as usize],
                    );
                } else {
                    self.add_vlb(ms, newi, new_vlb_col, vlb.coef, vlb.constant);
                }
            });

            // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 572-584 at f7b87ae01
            if keep {
                let i_hi = i as HighsInt;
                if detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, newi)) {
                    for val in 0..2 {
                        let idx = (2 * i_hi + val) as usize;
                        let mut entries: Vec<(HighsInt, HighsImplications_Implication)> =
                            Vec::new();
                        old_implication_tree[idx].for_each(&mut |k, v| entries.push((*k, *v)));
                        for (impl_col, impl_) in entries {
                            let new_impl_col = orig2reducedcol[impl_col as usize];
                            if new_impl_col == -1
                                || !ms
                                    .mipdata_
                                    .post_solve_stack
                                    .is_col_linearly_transformable(new_impl_col)
                            {
                                continue;
                            }
                            self.add_implication(newi, val == 1, new_impl_col, impl_);
                        }
                    }
                }
            }

            // todo also add old implications once implications can be added
            // incrementally for now we discard the old implications as they might be
            // weaker then newly computed ones and adding them would block computation
            // of new implications
        }
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::buildFrom (lines 561-581)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::buildFrom (lines 630-654)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::buildFrom (lines 588-612)
    /// C++: `void buildFrom(const HighsImplications & init)`
    pub fn build_from(&mut self, ms: &mut HighsMipSolver, init: &HighsImplications) {
        // todo check if this should be done
        let numcol = ms.num_col();

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 603-610 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;

        for (i, (init_vubs, init_vlbs)) in init
            .vubs
            .iter()
            .zip(init.vlbs.iter())
            .take(numcol as usize)
            .enumerate()
        {
            let i = i as HighsInt;
            init_vubs.for_each(&mut |&vub_col, vub| {
                let vub = *vub;
                if !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, vub_col)) {
                    return;
                }
                if vub.origin >= 0 {
                    self.add_vub_o(ms, i, vub_col, vub.coef, vub.constant, vub.origin);
                } else {
                    self.add_vub(ms, i, vub_col, vub.coef, vub.constant);
                }
            });

            init_vlbs.for_each(&mut |&vlb_col, vlb| {
                let vlb = *vlb;
                if !detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, vlb_col)) {
                    return;
                }
                if vlb.origin >= 0 {
                    self.add_vlb_o(ms, i, vlb_col, vlb.coef, vlb.constant, vlb.origin);
                } else {
                    self.add_vlb(ms, i, vlb_col, vlb.coef, vlb.constant);
                }
            });

            // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 603-610 at f7b87ae01
            if keep {
                if detach!(ms.mipdata_.domains[0], |gd| gd.is_binary(ms, i)) {
                    for val in 0..2 {
                        let idx = (2 * i + val) as usize;
                        let mut entries: Vec<(HighsInt, HighsImplications_Implication)> =
                            Vec::new();
                        init.implication_tree[idx].for_each(&mut |k, v| entries.push((*k, *v)));
                        for (impl_col, implic) in entries {
                            self.add_implication(i, val == 1, impl_col, implic);
                        }
                    }
                }
            }

            // todo also add old implications once implications can be added
            // incrementally for now we discard the old implications as they might be
            // weaker then newly computed ones and adding them would block computation
            // of new implications
        }
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::separateImpliedBounds (lines 583-752)
    /// C++: `void separateImpliedBounds(const HighsLpRelaxation & lpRelaxation, const std::vector<double> & sol, HighsCutPool & cutpool, double feastol, HighsDomain & globaldom, bool thread_safe)`
    ///
    /// PORT NOTE: signature changed because the C++ `globaldom` is `mipworker_.getGlobalDomain()`
    /// (HighsSeparation.cpp:94-97), which in the serial path is `mipdata_->getDomain()` itself:
    /// the same object `runProbing` changes through `mipsolver`. A detached `&mut HighsDomain`
    /// cannot alias it, so the domain is named by its index in `ms.mipdata_.domains` (the
    /// worker's `globaldom_` index, M2) and every read goes through `ms`, as the C++ aliasing does.
    pub fn separate_implied_bounds(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        sol: &[f64],
        cutpool: &mut HighsCutPool,
        feastol: f64,
        globaldom: usize,
        thread_safe: bool,
    ) {
        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 632, 672-731 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        let mut inds: [HighsInt; 2] = [0; 2];
        let mut vals: [f64; 2] = [0.0; 2];

        let mut numboundchgs: HighsInt = 0;

        // first do probing on all candidates that have not been probed yet
        if !ms.mipdata_.cliquetable.is_full() && !thread_safe {
            let oldnumqueries = ms.mipdata_.cliquetable.num_neighbourhood_queries;
            let old_num_entries: HighsInt = ms.mipdata_.cliquetable.get_num_entries();

            for &(col, _) in lp_relaxation.get_fractional_integers() {
                if ms.mipdata_.domains[globaldom].col_lower_[col as usize] != 0.0
                    || ms.mipdata_.domains[globaldom].col_upper_[col as usize] != 1.0
                    || (if keep {
                        self.probed_before(col, false) && self.probed_before(col, true)
                    } else {
                        self.implications_cached(col, false) && self.implications_cached(col, true)
                    })
                {
                    continue;
                }

                ms.profiling_.start(kMipClockProbingImplications, false);
                // PORT NOTE: the caller detached `cutpool`; probing propagates the global domain,
                // which replays that pool's events, so the pool is lent back for the call.
                let probing_result = crate::lend_in!(
                    ms.mipdata_.cutpools,
                    cutpool,
                    self.run_probing(ms, col, &mut numboundchgs)
                );
                ms.profiling_.stop(kMipClockProbingImplications);
                if probing_result && ms.mipdata_.domains[globaldom].infeasible() {
                    return;
                }

                if ms.mipdata_.cliquetable.is_full() {
                    break;
                }
            }

            let num_new_entries: HighsInt =
                ms.mipdata_.cliquetable.get_num_entries() - old_num_entries;

            self.next_cleanup_call -= cpp_max(0, num_new_entries);

            if self.next_cleanup_call < 0 {
                if !ms.parallel_lock_active() {
                    lend!(
                        ms.mipdata_.implications,
                        self,
                        crate::lend_in!(
                            ms.mipdata_.cutpools,
                            cutpool,
                            detach!(ms.mipdata_.cliquetable, |ct| detach!(
                                ms.mipdata_.domains[globaldom],
                                |gd| ct.run_clique_merging_2(ms, gd)
                            ))
                        )
                    );
                }

                self.next_cleanup_call = cpp_min(
                    ms.mipdata_.num_clique_entries_after_first_presolve,
                    ms.mipdata_.cliquetable.get_num_entries(),
                );
            }

            if !ms.parallel_lock_active() {
                ms.mipdata_.cliquetable.num_neighbourhood_queries = oldnumqueries;
            }
        }

        for &(col, _) in lp_relaxation.get_fractional_integers() {
            // skip non binary variables
            if ms.mipdata_.domains[globaldom].col_lower_[col as usize] != 0.0
                || ms.mipdata_.domains[globaldom].col_upper_[col as usize] != 1.0
            {
                continue;
            }

            // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 686-731 at f7b87ae01
            if keep {
                for val in 0..2 {
                    let mut entries: Vec<(HighsInt, HighsImplications_Implication)> = Vec::new();
                    self.implication_tree[(2 * col + val) as usize]
                        .for_each(&mut |k, v| entries.push((*k, *v)));
                    for (implcol, implic) in entries.into_iter() {
                        if val == 1 {
                            if implic.ub + feastol
                                < ms.mipdata_.domains[globaldom].col_upper_[implcol as usize]
                            {
                                vals[0] = 1.0;
                                inds[0] = implcol;
                                vals[1] = ms.mipdata_.domains[globaldom].col_upper_
                                    [implcol as usize]
                                    - implic.ub;
                                inds[1] = col;
                                let rhs =
                                    ms.mipdata_.domains[globaldom].col_upper_[implcol as usize];
                                Self::try_add_implied_cut(
                                    ms, cutpool, sol, feastol, &mut inds, &mut vals, rhs, implcol,
                                );
                            }
                            if implic.lb - feastol
                                > ms.mipdata_.domains[globaldom].col_lower_[implcol as usize]
                            {
                                vals[0] = -1.0;
                                inds[0] = implcol;
                                vals[1] = implic.lb
                                    - ms.mipdata_.domains[globaldom].col_lower_[implcol as usize];
                                inds[1] = col;
                                let rhs =
                                    -ms.mipdata_.domains[globaldom].col_lower_[implcol as usize];
                                Self::try_add_implied_cut(
                                    ms, cutpool, sol, feastol, &mut inds, &mut vals, rhs, implcol,
                                );
                            }
                        } else {
                            if implic.ub + feastol
                                < ms.mipdata_.domains[globaldom].col_upper_[implcol as usize]
                            {
                                vals[0] = 1.0;
                                inds[0] = implcol;
                                vals[1] = implic.ub
                                    - ms.mipdata_.domains[globaldom].col_upper_[implcol as usize];
                                inds[1] = col;
                                let rhs = implic.ub;
                                Self::try_add_implied_cut(
                                    ms, cutpool, sol, feastol, &mut inds, &mut vals, rhs, implcol,
                                );
                            }
                            if implic.lb - feastol
                                > ms.mipdata_.domains[globaldom].col_lower_[implcol as usize]
                            {
                                vals[0] = -1.0;
                                inds[0] = implcol;
                                vals[1] = ms.mipdata_.domains[globaldom].col_lower_
                                    [implcol as usize]
                                    - implic.lb;
                                inds[1] = col;
                                let rhs = -implic.lb;
                                Self::try_add_implied_cut(
                                    ms, cutpool, sol, feastol, &mut inds, &mut vals, rhs, implcol,
                                );
                            }
                        }
                    }
                }
                continue;
            }

            let mut infeas = false;
            if self.implications_cached(col, true) {
                let implics = crate::lend_in!(
                    ms.mipdata_.cutpools,
                    cutpool,
                    self.get_implications(ms, col, true, &mut infeas).to_vec()
                );
                if ms.mipdata_.domains[globaldom].infeasible() {
                    return;
                }

                if infeas {
                    vals[0] = 1.0;
                    inds[0] = col;
                    cutpool.add_cut(ms, &mut inds, &mut vals, 1, 0.0, false, true, false, false);
                    continue;
                }

                for chg in implics.iter() {
                    let ic = chg.column;
                    let ib = chg.boundval;
                    let rhs = if chg.boundtype == HighsBoundType::kUpper {
                        if ib + feastol >= ms.mipdata_.domains[globaldom].col_upper_[ic as usize] {
                            continue;
                        }

                        vals[0] = 1.0;
                        inds[0] = ic;
                        vals[1] = ms.mipdata_.domains[globaldom].col_upper_[ic as usize] - ib;
                        inds[1] = col;
                        ms.mipdata_.domains[globaldom].col_upper_[ic as usize]
                    } else {
                        if ib - feastol <= ms.mipdata_.domains[globaldom].col_lower_[ic as usize] {
                            continue;
                        }

                        vals[0] = -1.0;
                        inds[0] = ic;
                        vals[1] = ib - ms.mipdata_.domains[globaldom].col_lower_[ic as usize];
                        inds[1] = col;
                        -ms.mipdata_.domains[globaldom].col_lower_[ic as usize]
                    };

                    let viol =
                        sol[inds[0] as usize] * vals[0] + sol[inds[1] as usize] * vals[1] - rhs;

                    if viol > feastol {
                        // printf("added implied bound cut to pool\n");
                        let integral = Self::implied_cut_is_integral(ms, ic, &vals, rhs);
                        cutpool.add_cut(
                            ms, &mut inds, &mut vals, 2, rhs, integral, false, false, false,
                        );
                    }
                }
            }

            if self.implications_cached(col, false) {
                let implics = crate::lend_in!(
                    ms.mipdata_.cutpools,
                    cutpool,
                    self.get_implications(ms, col, false, &mut infeas).to_vec()
                );
                if ms.mipdata_.domains[globaldom].infeasible() {
                    return;
                }

                if infeas {
                    vals[0] = -1.0;
                    inds[0] = col;
                    cutpool.add_cut(ms, &mut inds, &mut vals, 1, -1.0, false, true, false, false);
                    continue;
                }

                for chg in implics.iter() {
                    let ic = chg.column;
                    let ib = chg.boundval;
                    let rhs = if chg.boundtype == HighsBoundType::kUpper {
                        if ib + feastol >= ms.mipdata_.domains[globaldom].col_upper_[ic as usize] {
                            continue;
                        }

                        vals[0] = 1.0;
                        inds[0] = ic;
                        vals[1] = ib - ms.mipdata_.domains[globaldom].col_upper_[ic as usize];
                        inds[1] = col;
                        ib
                    } else {
                        if ib - feastol <= ms.mipdata_.domains[globaldom].col_lower_[ic as usize] {
                            continue;
                        }

                        vals[0] = -1.0;
                        inds[0] = ic;
                        vals[1] = ms.mipdata_.domains[globaldom].col_lower_[ic as usize] - ib;
                        inds[1] = col;
                        -ib
                    };

                    let viol =
                        sol[inds[0] as usize] * vals[0] + sol[inds[1] as usize] * vals[1] - rhs;

                    if viol > feastol {
                        // printf("added implied bound cut to pool\n");
                        let integral = Self::implied_cut_is_integral(ms, ic, &vals, rhs);
                        cutpool.add_cut(
                            ms, &mut inds, &mut vals, 2, rhs, integral, false, false, false,
                        );
                    }
                }
            }
        }
    }

    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::separateImpliedBounds (lines 672-678)
    fn try_add_implied_cut(
        ms: &mut HighsMipSolver,
        cutpool: &mut HighsCutPool,
        sol: &[f64],
        feastol: f64,
        inds: &mut [HighsInt; 2],
        vals: &mut [f64; 2],
        rhs: f64,
        implcol: HighsInt,
    ) {
        let viol = sol[inds[0] as usize] * vals[0] + sol[inds[1] as usize] * vals[1] - rhs;
        if viol > feastol {
            let integral = Self::implied_cut_is_integral(ms, implcol, vals, rhs);
            cutpool.add_cut(ms, inds, vals, 2, rhs, integral, false, false, false);
        }
    }

    // PORT NOTE (correctness): HiGHS flags every
    // implied-bound cut on a non-continuous column as integral (`!isColContinuous(col)`). That flag
    // tells `HighsTransformedLp` the cut's slack is an integer variable, which is only true when
    // the binary's coefficient and the rhs are integers too. Implied bounds on implicit integers
    // can be fractional (x >= 1.93 * y), and an integral-flagged row with rhs -1.93 makes the
    // tableau CMIR cut off the optimum: 1.16.0 (#3320 rounding) exposes it on int_production.
    fn implied_cut_is_integral(
        ms: &HighsMipSolver,
        implcol: HighsInt,
        vals: &[f64; 2],
        rhs: f64,
    ) -> bool {
        let eps = ms.mipdata_.epsilon;
        let is_int = |v: f64| (v - v.round()).abs() <= eps;
        !ms.is_col_continuous(implcol) && is_int(vals[0]) && is_int(vals[1]) && is_int(rhs)
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::cleanupVarbounds (lines 754-793)
    /// C++: `void cleanupVarbounds(HighsInt col)`
    pub fn cleanup_varbounds(&mut self, ms: &mut HighsMipSolver, col: HighsInt) {
        let ub = ms.mipdata_.domains[0].col_upper_[col as usize];
        let lb = ms.mipdata_.domains[0].col_lower_[col as usize];

        if ub == lb {
            let mut num_vubs: i64 = 0;
            self.vubs[col as usize].for_each(&mut |_, _| num_vubs += 1);
            let mut num_vlbs: i64 = 0;
            self.vlbs[col as usize].for_each(&mut |_, _| num_vlbs += 1);
            self.num_var_bounds -= num_vubs + num_vlbs;
            self.vlbs[col as usize].clear();
            self.vubs[col as usize].clear();
            return;
        }

        let mut del_vbds: Vec<HighsInt> = Vec::new();
        // PORT NOTE (perf): the C++ lambdas mutate each `VarBound` in
        // place. The port does the same with `for_each_mut` while no entry needs the bound change
        // of `cleanupVub`/`cleanupVlb` (only that branch reads mutable solver state, and
        // `cleanup_vub` lends `self` back for it). From the first entry that needs it, the rest
        // of the tree runs on a snapshot (reused `vb_snapshot_` buffer) with the full
        // `cleanup_vub`/`cleanup_vlb` and a `find_mut` write-back, so entries are processed in
        // the same order with the same state as before.
        let mut snapshot = std::mem::take(&mut self.vb_snapshot_);

        let done = cleanup_vbds_in_place(
            &mut self.vubs[col as usize],
            ms,
            col,
            ub,
            true,
            &mut del_vbds,
        );
        if let Some(done) = done {
            collect_var_bounds_into(&self.vubs[col as usize], &mut snapshot);
            for &(vub_col, mut vub) in snapshot[done..].iter() {
                let mut redundant = false;
                let mut infeasible = false;
                self.cleanup_vub(
                    ms,
                    col,
                    vub_col,
                    &mut vub,
                    ub,
                    &mut redundant,
                    &mut infeasible,
                    true,
                );
                write_back_var_bound(&mut self.vubs[col as usize], vub_col, vub);
                if redundant {
                    del_vbds.push(vub_col);
                }
                // PORT NOTE: the C++ lambda `if (infeasible) return;` only ends this entry.
            }
        }

        for vub_col in del_vbds.iter() {
            self.vubs[col as usize].erase(vub_col);
        }
        self.num_var_bounds -= del_vbds.len() as i64;
        del_vbds.clear();

        let done = cleanup_vbds_in_place(
            &mut self.vlbs[col as usize],
            ms,
            col,
            lb,
            false,
            &mut del_vbds,
        );
        if let Some(done) = done {
            collect_var_bounds_into(&self.vlbs[col as usize], &mut snapshot);
            for &(vlb_col, mut vlb) in snapshot[done..].iter() {
                let mut redundant = false;
                let mut infeasible = false;
                self.cleanup_vlb(
                    ms,
                    col,
                    vlb_col,
                    &mut vlb,
                    lb,
                    &mut redundant,
                    &mut infeasible,
                    true,
                );
                write_back_var_bound(&mut self.vlbs[col as usize], vlb_col, vlb);
                if redundant {
                    del_vbds.push(vlb_col);
                }
                // PORT NOTE: the C++ lambda `if (infeasible) return;` only ends this entry.
            }
        }

        for vlb_col in del_vbds.iter() {
            self.vlbs[col as usize].erase(vlb_col);
        }
        self.num_var_bounds -= del_vbds.len() as i64;
        self.vb_snapshot_ = snapshot;
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::cleanupVlb (lines 795-834)
    /// C++: `void cleanupVlb(HighsInt col, HighsInt vlbCol, HighsImplications::VarBound & vlb, double lb, bool & redundant, bool & infeasible, bool allowBoundChanges = true) const`
    /// C++ default arguments: `allow_bound_changes = true`
    // PORT NOTE: `&mut self` (C++ `const`): the bound change lends `self` back to
    // `ms.mipdata_.implications` for `applyImplications` (re-entrant objects).
    pub fn cleanup_vlb(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vlb_col: HighsInt,
        vlb: &mut HighsImplications_VarBound,
        lb: f64,
        redundant: &mut bool,
        infeasible: &mut bool,
        allow_bound_changes: bool,
    ) {
        // initialize
        *redundant = false;
        *infeasible = false;

        // return if there is no variable bound
        if vlb_col == -1 {
            return;
        }

        let feastol = ms.mipdata_.feastol;
        let epsilon = ms.mipdata_.epsilon;

        // check variable lower bound
        ms.mipdata_
            .debug_solution
            .check_vlb(col, vlb_col, vlb.coef, vlb.constant);

        // HiGHS #3264 (mip_aggregate_implications): HighsImplications.cpp lines 977-1007 at
        // f2ae2b320 go through redundantVlb/tightenVlb, i.e. the inf-safe `VarBound` values.
        let inf_safe = ms.options_mip_.mip_aggregate_implications;
        let maxlb = HighsCDouble::from(var_bound_max_value_x(vlb, inf_safe));
        let minlb = HighsCDouble::from(var_bound_min_value_x(vlb, inf_safe));

        if maxlb <= lb + feastol {
            // variable bound is redundant
            *redundant = true;
        } else if minlb < lb - epsilon {
            // coefficient can be tightened
            let newcoef = f64::from(HighsCDouble::from(lb) - maxlb);
            if vlb.coef < 0.0 {
                vlb.coef = newcoef;
            } else {
                vlb.constant = lb;
                vlb.coef = -newcoef;
            }
            // check tightened variable lower bound
            ms.mipdata_
                .debug_solution
                .check_vlb(col, vlb_col, vlb.coef, vlb.constant);
        } else if allow_bound_changes && minlb > lb + epsilon {
            let boundval = f64::from(minlb);
            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                    ms,
                    HighsDomainChange {
                        boundval,
                        column: col,
                        boundtype: HighsBoundType::kLower,
                    },
                    HighsDomain_Reason::unspecified(),
                ))
            );
            *infeasible = ms.mipdata_.domains[0].infeasible();
        }
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::cleanupVub (lines 836-875)
    /// C++: `void cleanupVub(HighsInt col, HighsInt vubCol, HighsImplications::VarBound & vub, double ub, bool & redundant, bool & infeasible, bool allowBoundChanges = true) const`
    /// C++ default arguments: `allow_bound_changes = true`
    // PORT NOTE: `&mut self` (C++ `const`): the bound change lends `self` back to
    // `ms.mipdata_.implications` for `applyImplications` (re-entrant objects).
    pub fn cleanup_vub(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        vub_col: HighsInt,
        vub: &mut HighsImplications_VarBound,
        ub: f64,
        redundant: &mut bool,
        infeasible: &mut bool,
        allow_bound_changes: bool,
    ) {
        // initialize
        *redundant = false;
        *infeasible = false;

        // return if there is no variable bound
        if vub_col == -1 {
            return;
        }

        let feastol = ms.mipdata_.feastol;
        let epsilon = ms.mipdata_.epsilon;

        // check variable upper bound
        ms.mipdata_
            .debug_solution
            .check_vub(col, vub_col, vub.coef, vub.constant);

        // HiGHS #3264 (mip_aggregate_implications): HighsImplications.cpp lines 1009-1039 at
        // f2ae2b320 go through redundantVub/tightenVub, i.e. the inf-safe `VarBound` values.
        let inf_safe = ms.options_mip_.mip_aggregate_implications;
        let maxub = HighsCDouble::from(var_bound_max_value_x(vub, inf_safe));
        let minub = HighsCDouble::from(var_bound_min_value_x(vub, inf_safe));

        if minub >= ub - feastol {
            // variable bound is redundant
            *redundant = true;
        } else if maxub > ub + epsilon {
            // coefficient can be tightened
            let newcoef = f64::from(HighsCDouble::from(ub) - minub);
            if vub.coef > 0.0 {
                vub.coef = newcoef;
            } else {
                vub.constant = ub;
                vub.coef = -newcoef;
            }
            // check tightened variable upper bound
            ms.mipdata_
                .debug_solution
                .check_vub(col, vub_col, vub.coef, vub.constant);
        } else if allow_bound_changes && maxub < ub - epsilon {
            let boundval = f64::from(maxub);
            lend!(
                ms.mipdata_.implications,
                self,
                detach!(ms.mipdata_.domains[0], |gd| gd.change_bound(
                    ms,
                    HighsDomainChange {
                        boundval,
                        column: col,
                        boundtype: HighsBoundType::kUpper,
                    },
                    HighsDomain_Reason::unspecified(),
                ))
            );
            *infeasible = ms.mipdata_.domains[0].infeasible();
        }
    }

    /// HiGHS: HighsImplications.cpp::HighsImplications::applyImplications (lines 877-916)
    /// C++: `void applyImplications(HighsDomain & domain, HighsInt col, HighsInt val)`
    ///
    /// PORT NOTE: signature changed because the body needs `ms` (`variableType`, `feastol`,
    /// `changeBound`), and because the call is re-entrant: `HighsDomain::changeBound`
    /// (HighsDomain.cpp:2100-2101) calls it, and its own `domain.changeBound` can call it again.
    /// It is therefore an associated function that reads `ms.mipdata_.implications` afresh at
    /// every step instead of holding a detached `&mut self`. The implication list is not
    /// written during the loop (C++ iterates it by reference), so indexing it by position is
    /// the same walk. Callers must not hold `ms.mipdata_.implications` detached (re-entrant
    /// objects).
    pub fn apply_implications(
        ms: &mut HighsMipSolver,
        domain: &mut HighsDomain,
        col: HighsInt,
        val: HighsInt,
    ) {
        debug_assert!(domain.is_fixed(col));

        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 855-886 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        if keep {
            let mut entries: Vec<(HighsInt, HighsImplications_Implication)> = Vec::new();
            ms.mipdata_.implications.implication_tree[(2 * col + val) as usize]
                .for_each(&mut |k, v| entries.push((*k, *v)));
            for (implcol, implic) in entries.into_iter() {
                if Self::check_implication_x1b(ms, domain, col, val, implcol, implic) {
                    break;
                }
            }
            return;
        }

        let loc = (2 * col + val) as usize;
        if ms.mipdata_.implications.implications[loc].computed {
            let mut k = 0usize;
            while k < ms.mipdata_.implications.implications[loc].implics.len() {
                let domchg = ms.mipdata_.implications.implications[loc].implics[k].clone();
                if Self::check_implication(ms, domain, &domchg, col, val) {
                    break;
                }
                k += 1;
            }
        }
    }

    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.cpp::HighsImplications::applyImplications (lines 857-883)
    fn check_implication_x1b(
        ms: &mut HighsMipSolver,
        domain: &mut HighsDomain,
        col: HighsInt,
        val: HighsInt,
        implcol: HighsInt,
        implic: HighsImplications_Implication,
    ) -> bool {
        debug_assert!(!domain.infeasible());
        if domain.is_fixed(implcol) {
            return false;
        }
        let isint = domain.variable_type(ms, implcol) != HighsVarType::kContinuous;
        // Directly change bounds on all integer columns. Only change continuous
        // columns that fix the column, as changing their domains risks
        // suppressing further bound changes found in propagation, e.g.,
        // change [0, 100] -> [0, 50], propagation could tighten to [0, 48], but
        // such a tightening would not be applied due to min boundRange improvement.
        let dcol = implcol as usize;
        if (!isint && implic.lb > domain.col_upper_[dcol] - domain.feastol(ms))
            || (isint && implic.lb > domain.col_lower_[dcol] + domain.feastol(ms))
        {
            domain.change_bound(
                ms,
                HighsDomainChange {
                    boundval: implic.lb,
                    column: implcol,
                    boundtype: HighsBoundType::kLower,
                },
                HighsDomain_Reason::clique_table(col, val),
            );
            if domain.infeasible() {
                return true;
            }
        }
        if (!isint && implic.ub < domain.col_lower_[dcol] + domain.feastol(ms))
            || (isint && implic.ub < domain.col_upper_[dcol] - domain.feastol(ms))
        {
            domain.change_bound(
                ms,
                HighsDomainChange {
                    boundval: implic.ub,
                    column: implcol,
                    boundtype: HighsBoundType::kUpper,
                },
                HighsDomain_Reason::clique_table(col, val),
            );
        }
        domain.infeasible()
    }

    /// PORT NOTE (lambda `checkImplication`, HighsImplications.cpp:881-909).
    fn check_implication(
        ms: &mut HighsMipSolver,
        domain: &mut HighsDomain,
        domchg: &HighsDomainChange,
        col: HighsInt,
        val: HighsInt,
    ) -> bool {
        debug_assert!(!domain.infeasible());
        if domain.is_fixed(domchg.column) {
            return false;
        }
        let isint = domain.variable_type(ms, domchg.column) != HighsVarType::kContinuous;
        // Directly change bounds on all integer columns. Only change continuous
        // columns that fix the column, as changing their domains risks
        // suppressing further bound changes found in propagation, e.g.,
        // change [0, 100] -> [0, 50], propagation could tighten to [0, 48], but
        // such a tightening would not be applied due to min boundRange improvement.
        let dcol = domchg.column as usize;
        if domchg.boundtype == HighsBoundType::kLower {
            if (!isint && domchg.boundval > domain.col_upper_[dcol] - domain.feastol(ms))
                || (isint && domchg.boundval > domain.col_lower_[dcol] + domain.feastol(ms))
            {
                domain.change_bound(
                    ms,
                    domchg.clone(),
                    HighsDomain_Reason::clique_table(col, val),
                );
            }
        } else if (!isint && domchg.boundval < domain.col_lower_[dcol] + domain.feastol(ms))
            || (isint && domchg.boundval < domain.col_upper_[dcol] - domain.feastol(ms))
        {
            domain.change_bound(
                ms,
                domchg.clone(),
                HighsDomain_Reason::clique_table(col, val),
            );
        }
        domain.infeasible()
    }

    /// HiGHS: HighsImplications.h::HighsImplications::HighsImplications (lines 60-70)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::HighsImplications (lines 80-85)
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::HighsImplications (lines 104-110, with resize lines 137-148)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications::HighsImplications (lines 73-78, with resize lines 102-110)
    /// C++: `HighsImplications(const HighsMipSolver & mipsolver)`
    pub fn new(mipsolver: &HighsMipSolver) -> Self {
        let numcol = mipsolver.num_col();
        let numcol_us = numcol as usize;
        let mut this = Self {
            next_cleanup_call: 0,
            implications: Vec::new(),
            num_implications: 0,
            num_var_bounds: 0,
            max_var_bounds: 0,
            vubs: Vec::new(),
            vlbs: Vec::new(),
            row_to_var_bounds: Vec::new(),
            substitutions: Vec::new(),
            colsubstituted: Vec::new(),
            store_lifting_opportunity: None,
            implications_down: Vec::new(),
            implications_up: Vec::new(),
            dual_fix_probing_bin_inds_: Vec::new(),
            dual_fix_probing_bin_flags_: Vec::new(),
            implication_tree: Vec::new(),
            reverse_implications: Vec::new(),
            has_probed: Vec::new(),
            keep_x1b: false,
            vb_snapshot_: Vec::new(),
        };
        this.implications
            .resize_with(2 * numcol_us, HighsImplications_Implics::default);
        this.colsubstituted.resize(numcol_us, 0);
        this.vubs.resize_with(numcol_us, HighsHashTree::default);
        this.vlbs.resize_with(numcol_us, HighsHashTree::default);
        this.resize_rows(mipsolver.num_row());
        this.next_cleanup_call = mipsolver.num_nonzero();
        this.num_implications = 0;
        this.num_var_bounds = 0;
        this.max_var_bounds = Self::calc_max_var_bounds(numcol);
        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.h lines 146-147 at 8a5f9c1d8
        if mipsolver.options_mip_.mip_dual_fix_probing {
            this.dual_fix_probing_bin_inds_.reserve(numcol_us);
            this.dual_fix_probing_bin_flags_ =
                vec![HighsImplications_TentativeFixing::default(); numcol_us];
        }
        // HiGHS #3209 (mip_keep_implications): HighsImplications.h lines 102-105 at f7b87ae01
        if mipsolver.options_mip_.mip_keep_implications {
            this.keep_x1b = true;
            this.implication_tree
                .resize_with(2 * numcol_us, HighsHashTree::default);
            this.has_probed.resize(2 * numcol_us, false);
            this.reverse_implications
                .resize_with(numcol_us, HighsHashTree::default);
        }
        this
    }

    /// HiGHS: HighsImplications.h::HighsImplications::reset (lines 75-95)
    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::reset (lines 90-109)
    /// HiGHS latest (#3244, 8a5f9c1d8): HighsImplications.h::HighsImplications::reset (lines 115-133, with resize lines 137-148)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications::reset (lines 83-100, with resize lines 102-110)
    /// C++: `void reset()`
    ///
    /// PORT NOTE: as in C++, `substitutions` is not cleared by reset.
    pub fn reset(&mut self, ms: &mut HighsMipSolver) {
        self.colsubstituted.clear();
        self.colsubstituted.shrink_to_fit();
        self.implications.clear();
        self.implications.shrink_to_fit();

        // HiGHS #3209 (mip_keep_implications): HighsImplications.h lines 88-91 at f7b87ae01
        let keep = ms.options_mip_.mip_keep_implications;
        if keep {
            self.implication_tree.clear();
            self.implication_tree.shrink_to_fit();
            self.has_probed.clear();
            self.has_probed.shrink_to_fit();
            self.reverse_implications.clear();
            self.reverse_implications.shrink_to_fit();
        }

        let numcol = ms.num_col();
        let numcol_us = numcol as usize;
        self.implications
            .resize_with(2 * numcol_us, HighsImplications_Implics::default);
        self.colsubstituted.resize(numcol_us, 0);
        self.num_implications = 0;
        self.vubs.clear();
        self.vubs.shrink_to_fit();
        self.vubs.resize_with(numcol_us, HighsHashTree::default);
        self.vlbs.clear();
        self.vlbs.shrink_to_fit();
        self.vlbs.resize_with(numcol_us, HighsHashTree::default);
        // HiGHS #3244 (mip_dual_fix_probing): HighsImplications.h lines 129, 146-147 at 8a5f9c1d8
        if ms.options_mip_.mip_dual_fix_probing {
            self.dual_fix_probing_bin_inds_.clear();
            self.dual_fix_probing_bin_inds_.reserve(numcol_us);
            self.dual_fix_probing_bin_flags_.clear();
            self.dual_fix_probing_bin_flags_
                .resize(numcol_us, HighsImplications_TentativeFixing::default());
        }
        // HiGHS #3209 (mip_keep_implications): HighsImplications.h lines 102-105 at f7b87ae01
        if keep {
            self.implication_tree
                .resize_with(2 * numcol_us, HighsHashTree::default);
            self.has_probed.resize(2 * numcol_us, false);
            self.reverse_implications
                .resize_with(numcol_us, HighsHashTree::default);
        }
        self.row_to_var_bounds.clear();
        self.row_to_var_bounds.shrink_to_fit();
        self.resize_rows(ms.num_row());
        self.num_var_bounds = 0;
        self.max_var_bounds = Self::calc_max_var_bounds(numcol);

        self.next_cleanup_call = ms.num_nonzero();
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::resize (lines 111-120)
    pub fn resize_rows(&mut self, nrows: HighsInt) {
        self.row_to_var_bounds
            .resize_with(nrows as usize, HighsHashTree::default);
    }

    /// HiGHS 1.16 (#3218): HighsImplications.h::HighsImplications::resize(ncols, nrows)
    /// (1.16.0 lines 137-148), called by `HPresolve::addToMatrix` (rows) under
    /// `mip_presolve_row_append`. The vectors that exist only under their own HiGHS backport/HiGHS backport
    /// option in this port (`implication_tree`/`has_probed`/`reverse_implications`,
    /// `dual_fix_probing_*`) are resized only when that option is on, as in `reset`.
    pub fn resize_x1b(&mut self, ncols: HighsInt, nrows: HighsInt, dual_fix_probing: bool) {
        let n = ncols as usize;
        self.implications
            .resize_with(2 * n, HighsImplications_Implics::default);
        if self.keep_x1b {
            self.implication_tree
                .resize_with(2 * n, HighsHashTree::default);
            self.has_probed.resize(2 * n, false);
            self.reverse_implications
                .resize_with(n, HighsHashTree::default);
        }
        self.colsubstituted.resize(n, 0);
        self.vubs.resize_with(n, HighsHashTree::default);
        self.vlbs.resize_with(n, HighsHashTree::default);
        self.resize_rows(nrows);
        self.max_var_bounds = Self::calc_max_var_bounds(ncols);
        if dual_fix_probing {
            self.dual_fix_probing_bin_inds_
                .reserve(n.saturating_sub(self.dual_fix_probing_bin_inds_.len()));
            self.dual_fix_probing_bin_flags_.clear();
            self.dual_fix_probing_bin_flags_
                .resize(n, HighsImplications_TentativeFixing::default());
        }
    }

    /// HiGHS 1.16 (#3250): HighsImplications.cpp::HighsImplications::reindexRows (1.16.0
    /// lines 635-655), called by `HPresolve::moveCutsToPool` under
    /// `mip_presolve_row_append`.
    pub fn reindex_rows(&mut self, nrows: HighsInt, new_row_index: &[HighsInt]) {
        let Self {
            row_to_var_bounds,
            vlbs,
            vubs,
            ..
        } = self;
        row_to_var_bounds.clear();
        row_to_var_bounds.resize_with(nrows as usize, HighsHashTree::default);

        for col in 0..vlbs.len() {
            let col_i = col as HighsInt;
            vlbs[col].for_each_mut(&mut |&vlb_col, vlb| {
                if vlb.origin >= 0 {
                    vlb.origin = new_row_index[vlb.origin as usize];
                    if vlb.origin >= 0 {
                        row_to_var_bounds[vlb.origin as usize]
                            .insert(HighsHashTableEntry::new(col_i, vlb_col));
                    }
                }
            });
            vubs[col].for_each_mut(&mut |&vub_col, vub| {
                if vub.origin >= 0 {
                    vub.origin = new_row_index[vub.origin as usize];
                    if vub.origin >= 0 {
                        row_to_var_bounds[vub.origin as usize]
                            .insert(HighsHashTableEntry::new(col_i, vub_col));
                    }
                }
            });
        }
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::getVlbs (lines 190-192)
    pub fn get_vlbs(&self, col: HighsInt) -> &HighsHashTree<HighsInt, HighsImplications_VarBound> {
        &self.vlbs[col as usize]
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::getVubs (lines 193-195)
    pub fn get_vubs(&self, col: HighsInt) -> &HighsHashTree<HighsInt, HighsImplications_VarBound> {
        &self.vubs[col as usize]
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.h::HighsImplications::getRowVarBounds (lines 199-201)
    pub fn get_row_var_bounds(&self, row: HighsInt) -> &HighsHashTree<HighsInt, HighsInt> {
        &self.row_to_var_bounds[row as usize]
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::rowModified (lines 415-423)
    pub fn row_modified(&mut self, row: HighsInt) {
        let mut entries: Vec<(HighsInt, HighsInt)> = Vec::new();
        self.row_to_var_bounds[row as usize].for_each(&mut |col, binary_col| {
            entries.push((*col, *binary_col));
        });
        for (col, binary_col) in entries {
            if let Some(vb) = self.vlbs[col as usize].find_mut(&binary_col) {
                if vb.origin == row {
                    vb.origin = -1;
                }
            }
            if let Some(vb) = self.vubs[col as usize].find_mut(&binary_col) {
                if vb.origin == row {
                    vb.origin = -1;
                }
            }
        }
        self.row_to_var_bounds[row as usize].clear();
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::redundantVlb (lines 815-817)
    pub fn redundant_vlb(
        &self,
        ms: &HighsMipSolver,
        vlb: &HighsImplications_VarBound,
        lb: f64,
    ) -> bool {
        var_bound_max_value_inf_safe(vlb) <= lb + ms.mipdata_.feastol
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::redundantVub (lines 819-821)
    pub fn redundant_vub(
        &self,
        ms: &HighsMipSolver,
        vub: &HighsImplications_VarBound,
        ub: f64,
    ) -> bool {
        var_bound_min_value_inf_safe(vub) >= ub - ms.mipdata_.feastol
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::tightenVlb (lines 823-836)
    pub fn tighten_vlb(
        &self,
        ms: &HighsMipSolver,
        vlb: &mut HighsImplications_VarBound,
        lb: f64,
    ) -> bool {
        if self.redundant_vlb(ms, vlb, lb)
            || var_bound_min_value_inf_safe(vlb) >= lb - ms.mipdata_.epsilon
        {
            return false;
        }
        let newcoef = f64::from(lb - HighsCDouble::from(var_bound_max_value_inf_safe(vlb)));
        if vlb.coef < 0.0 {
            vlb.coef = newcoef;
        } else {
            vlb.constant = lb;
            vlb.coef = -newcoef;
        }
        true
    }

    /// HiGHS latest (#3264, f2ae2b320): HighsImplications.cpp::HighsImplications::tightenVub (lines 838-851)
    pub fn tighten_vub(
        &self,
        ms: &HighsMipSolver,
        vub: &mut HighsImplications_VarBound,
        ub: f64,
    ) -> bool {
        if self.redundant_vub(ms, vub, ub)
            || var_bound_max_value_inf_safe(vub) <= ub + ms.mipdata_.epsilon
        {
            return false;
        }
        let newcoef = f64::from(ub - HighsCDouble::from(var_bound_min_value_inf_safe(vub)));
        if vub.coef > 0.0 {
            vub.coef = newcoef;
        } else {
            vub.constant = ub;
            vub.coef = -newcoef;
        }
        true
    }

    /// HiGHS: HighsImplications.h::HighsImplications::calcMaxVarBounds (lines 97-99)
    /// C++: `static int64_t calcMaxVarBounds(HighsInt numcol)`
    pub fn calc_max_var_bounds(numcol: HighsInt) -> i64 {
        5_000_000 + 10 * numcol as i64
    }

    /// HiGHS: HighsImplications.h::HighsImplications::getNumImplications (lines 101-103)
    /// C++: `HighsInt getNumImplications() const`
    pub fn get_num_implications(&self) -> HighsInt {
        self.num_implications as HighsInt
    }

    /// HiGHS: HighsImplications.h::HighsImplications::getImplications (lines 105-116)
    /// C++: `const std::vector<HighsDomainChange> & getImplications(HighsInt col, bool val, bool & infeasible)`
    pub fn get_implications(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: bool,
        infeasible: &mut bool,
    ) -> &[HighsDomainChange] {
        let loc = (2 * col + HighsInt::from(val)) as usize;
        if !self.implications[loc].computed {
            *infeasible = self.compute_implications(ms, col, val);
        } else {
            *infeasible = false;
        }

        debug_assert!(self.implications[loc].computed);

        &self.implications[loc].implics
    }

    /// HiGHS: HighsImplications.h::HighsImplications::implicationsCached (lines 118-121)
    /// C++: `bool implicationsCached(HighsInt col, bool val)`
    pub fn implications_cached(&mut self, col: HighsInt, val: bool) -> bool {
        let loc = (2 * col + HighsInt::from(val)) as usize;
        self.implications[loc].computed
    }

    /// HiGHS: HighsImplications.h::HighsImplications::tooManyVarBounds (lines 123-123)
    /// C++: `bool tooManyVarBounds() const`
    pub fn too_many_var_bounds(&self) -> bool {
        self.num_var_bounds >= self.max_var_bounds
    }

    /// HiGHS: HighsImplications.h::HighsImplications::columnTransformed (lines 139-160)
    /// HiGHS latest (#3209, f7b87ae01): HighsImplications.h::HighsImplications::columnTransformed (lines 142-156)
    /// C++: `void columnTransformed(HighsInt col, double scale, double constant)`
    pub fn column_transformed(&mut self, col: HighsInt, scale: f64, constant: f64) {
        let ci = col as usize;

        // HiGHS #3209 (mip_keep_implications): HighsImplications.h lines 143-156 at f7b87ae01
        if self.keep_x1b {
            let mut bin_cols: Vec<HighsInt> = Vec::new();
            self.reverse_implications[ci].for_each(&mut |k, _| bin_cols.push(*k));
            for bin_col in bin_cols {
                for val in 0..2 {
                    let idx = (2 * bin_col + val) as usize;
                    if let Some(implic) = self.implication_tree[idx].find_mut(&col) {
                        if scale < 0.0 {
                            std::mem::swap(&mut implic.lb, &mut implic.ub);
                        }
                        implic.lb -= constant;
                        implic.lb /= scale;
                        implic.ub -= constant;
                        implic.ub /= scale;
                    }
                }
            }
        }

        // Update variable bounds affected by transformation
        if scale < 0.0 {
            std::mem::swap(&mut self.vubs[ci], &mut self.vlbs[ci]);
        }

        let transform = |vbd: &mut HighsImplications_VarBound| {
            vbd.constant -= constant;
            vbd.constant /= scale;
            vbd.coef /= scale;
        };

        self.vlbs[ci].for_each_mut(&mut |_, vbd| transform(vbd));
        self.vubs[ci].for_each_mut(&mut |_, vbd| transform(vbd));

        // Update substitutions affected by transformation
        for substitution in self.substitutions.iter_mut() {
            if substitution.substcol == col {
                substitution.offset -= constant;
                substitution.offset /= scale;
                substitution.scale /= scale;
            }
        }
    }
}

/// HiGHS: HighsImplications.cpp::HighsImplications::getBestVub, lambda `isVubBetter` (lines 168-185)
///
/// PORT NOTE: C++ indexes `lpSolution.col_dual[bestVub.first]` even when `bestVub.first == -1`
/// (undefined behaviour). That state needs an infinite (or NaN) `ub_dist`, `vub_nodes == 0` and
/// `min_vub_val` within feastol of `minbest_ub`. The C++ term is then `col_dual[-1] / 0.0`,
/// because the initial `bestVub.second.coef` is 0.0. For the usual heap layout `col_dual[-1]` is a
/// positive denormal (the allocator size word), which gives `x - inf > feastol` = false. A zero
/// gives NaN, which is also false. We return false without indexing with -1.
#[allow(clippy::too_many_arguments)]
fn is_vub_better(
    ub_dist: f64,
    vub_nodes: i64,
    min_vub_val: f64,
    vub_col: HighsInt,
    vub: &HighsImplications_VarBound,
    best: &(HighsInt, HighsImplications_VarBound),
    best_ub_dist: f64,
    bestvubnodes: i64,
    minbest_ub: f64,
    lp_solution: &HighsSolution,
    feastol: f64,
) -> bool {
    if ub_dist < best_ub_dist - feastol {
        return true;
    }
    if vub_nodes > bestvubnodes {
        return true;
    }
    if vub_nodes < bestvubnodes {
        return false;
    }
    if min_vub_val < minbest_ub - feastol {
        return true;
    }
    if K_USE_DUALS_FOR_BREAKING_TIES {
        if min_vub_val > minbest_ub + feastol {
            return false;
        }
        if best.0 < 0 {
            return false;
        }
        if lp_solution.col_dual[vub_col as usize] / vub.coef
            - lp_solution.col_dual[best.0 as usize] / best.1.coef
            > feastol
        {
            return true;
        }
    }

    false
}

/// HiGHS: HighsImplications.cpp::HighsImplications::getBestVlb, lambda `isVlbBetter` (lines 242-258)
///
/// PORT NOTE: see [`is_vub_better`] for the `best.0 == -1` tie-break.
#[allow(clippy::too_many_arguments)]
fn is_vlb_better(
    lb_dist: f64,
    vlb_nodes: i64,
    max_vlb_val: f64,
    vlb_col: HighsInt,
    vlb: &HighsImplications_VarBound,
    best: &(HighsInt, HighsImplications_VarBound),
    best_lb_dist: f64,
    bestvlbnodes: i64,
    maxbest_lb: f64,
    lp_solution: &HighsSolution,
    feastol: f64,
) -> bool {
    if lb_dist < best_lb_dist - feastol {
        return true;
    }
    if vlb_nodes > bestvlbnodes {
        return true;
    }
    if vlb_nodes < bestvlbnodes {
        return false;
    }
    if max_vlb_val > maxbest_lb + feastol {
        return true;
    }
    if K_USE_DUALS_FOR_BREAKING_TIES {
        if max_vlb_val < maxbest_lb - feastol {
            return false;
        }
        if best.0 < 0 {
            return false;
        }
        if lp_solution.col_dual[vlb_col as usize] / vlb.coef
            - lp_solution.col_dual[best.0 as usize] / best.1.coef
            < -feastol
        {
            return true;
        }
    }

    false
}

impl HighsImplications_VarBound {
    /// HiGHS: HighsImplications.h::HighsImplications::VarBound::minValue (lines 42-45)
    /// C++: `double minValue() const`
    pub fn min_value(&self) -> f64 {
        var_bound_min_value(self)
    }

    /// HiGHS: HighsImplications.h::HighsImplications::VarBound::maxValue (lines 46-49)
    /// C++: `double maxValue() const`
    pub fn max_value(&self) -> f64 {
        var_bound_max_value(self)
    }

    /// `minValue()` of the header in force: the 1.15.1 one, or with `inf_safe` the #3264 one
    /// (f2ae2b320, HighsImplications.h lines 56-61), which 1.16.0 uses everywhere.
    pub fn min_value_x(&self, inf_safe: bool) -> f64 {
        var_bound_min_value_x(self, inf_safe)
    }

    /// `maxValue()` of the header in force, see [`Self::min_value_x`].
    pub fn max_value_x(&self, inf_safe: bool) -> f64 {
        var_bound_max_value_x(self, inf_safe)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn x1b_impl_b_probing_reuses_cached() {
        // HiGHS #3209 (mip_keep_implications): HighsImplications.cpp lines 305-307 at f7b87ae01
        let mut ms = x1b_impl_a_model(true);
        let mut nred: HighsInt = 0;
        let first = detach!(ms.mipdata_.implications, |imp| imp
            .run_probing(&mut ms, 0, &mut nred));
        assert!(first);
        assert!(ms.mipdata_.implications.probed_before(0, true));
        assert!(ms.mipdata_.implications.probed_before(0, false));
        let n1 = ms.mipdata_.implications.get_num_implications();
        let second = detach!(ms.mipdata_.implications, |imp| imp
            .run_probing(&mut ms, 0, &mut nred));
        assert!(!second);
        assert_eq!(ms.mipdata_.implications.get_num_implications(), n1);
    }
    fn x1b_impl_a_model(keep: bool) -> HighsMipSolver {
        // same LP as x1_dfp_c_probe_model, with the HiGHS #3209 option set before construction
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 4;
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 10.0, 1.0];
        lp.col_cost_ = vec![0.0, 1.0, 0.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = vec![2.0, -kHighsInf, 1.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, 10.0, kHighsInf, 10.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 6, 8];
        lp.a_matrix_.index_ = vec![0, 1, 2, 0, 1, 3, 2, 3];
        lp.a_matrix_.value_ = vec![2.0, 8.0, 1.0, 1.0, 1.0, 1.0, 1.0, 6.0];

        let mut options = HighsOptions::default();
        options.mip_keep_implications = keep;
        options.mip_cont_prop_improve = 0.3; // pins HiGHS threshold
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.setup_domain_propagation();
        // rebuild asks the postsolve stack whether a column is linearly transformable
        ms.mipdata_.post_solve_stack.initialize_index_maps(4, 3);
        ms
    }

    #[test]
    fn x1b_impl_a_rebuild_keeps_implications() {
        for keep in [false, true] {
            let mut ms = x1b_impl_a_model(keep);
            let infeasible = detach!(ms.mipdata_.implications, |imp| imp
                .compute_implications(&mut ms, 0, true));
            assert!(!infeasible);
            assert!(ms.mipdata_.implications.get_num_implications() > 0);
            detach!(ms.mipdata_.implications, |imp| imp.rebuild(
                &mut ms,
                3,
                &[0, 1, 2],
                &[0, 1, 2, 3]
            ));
            let n = ms.mipdata_.implications.get_num_implications();
            if keep {
                assert!(n > 0);
            } else {
                assert_eq!(n, 0);
            }
        }
    }

    #[test]
    fn x1b_impl_a_column_transformed_updates() {
        let mut ms = x1b_impl_a_model(true);
        ms.mipdata_.implications.add_implication(
            0,
            true,
            1,
            HighsImplications_Implication { lb: 2.0, ub: 6.0 },
        );
        ms.mipdata_.implications.column_transformed(1, -2.0, 1.0);
        let entry = ms.mipdata_.implications.implication_tree[1]
            .find(&1)
            .copied()
            .expect("entry must exist");
        // swapped: lb = 6, ub = 2; then (x - 1) / -2
        assert_eq!(entry.lb, -2.5);
        assert_eq!(entry.ub, -0.5);
    }

    /// x0, x2 binary, x1 continuous in [0, 10]; rows
    /// r0: 2 x0 + x1 >= 2, r1: 8 x0 + x1 <= 10, r2: x0 + x2 >= 1, r3: x1 + 6 x2 <= 10.
    /// No row tightens x1 globally. Probing x0: down (x0 = 0) gives x2 >= 1 (binary, a
    /// clique) and through r3 x1 in [2, 4]; up (x0 = 1) gives x1 <= 2. The merge then
    /// tightens the global x1 <= 4.
    fn x1_dfp_c_probe_model(dual_fix_probing: bool) -> HighsMipSolver {
        x1_dfp_c_probe_model_keep(dual_fix_probing, false)
    }

    fn x1_dfp_c_probe_model_keep(dual_fix_probing: bool, keep: bool) -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 3;
        lp.num_row_ = 4;
        lp.col_lower_ = vec![0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 10.0, 1.0];
        lp.col_cost_ = vec![0.0, 1.0, 0.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
            HighsVarType::kInteger,
        ];
        lp.row_lower_ = vec![2.0, -kHighsInf, 1.0, -kHighsInf];
        lp.row_upper_ = vec![kHighsInf, 10.0, kHighsInf, 10.0];
        // column-wise: x0 in r0 (2), r1 (8), r2 (1); x1 in r0 (1), r1 (1), r3 (1);
        // x2 in r2 (1), r3 (6)
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 3;
        lp.a_matrix_.num_row_ = 4;
        lp.a_matrix_.start_ = vec![0, 3, 6, 8];
        lp.a_matrix_.index_ = vec![0, 1, 2, 0, 1, 3, 2, 3];
        lp.a_matrix_.value_ = vec![2.0, 8.0, 1.0, 1.0, 1.0, 1.0, 1.0, 6.0];

        let mut options = HighsOptions::default();
        options.mip_dual_fix_probing = dual_fix_probing;
        options.mip_keep_implications = keep;
        options.mip_cont_prop_improve = 0.3; // pins HiGHS threshold
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.new_data();
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.setup_domain_propagation();
        ms
    }

    fn x1_dfp_c_chg(
        column: HighsInt,
        boundtype: HighsBoundType,
        boundval: f64,
    ) -> (HighsInt, i32, f64) {
        (column, boundtype as i32, boundval)
    }

    fn x1_dfp_c_key(v: &[HighsDomainChange]) -> Vec<(HighsInt, i32, f64)> {
        v.iter()
            .map(|chg| x1_dfp_c_chg(chg.column, chg.boundtype, chg.boundval))
            .collect()
    }

    #[test]
    fn x1_dfp_c_merge_keeps_old_order_off() {
        // Option off: runProbing merges the 1.15.1 getImplications lists (sorted by
        // HighsDomainChange::operator<) and the new vectors stay empty.
        let mut ms = x1_dfp_c_probe_model(false);
        assert!(!ms.mipdata_.domains[0].get_dual_fix_probing_active());
        let mut num_reductions: HighsInt = 0;
        let probed = detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_reductions
        ));
        assert!(probed);
        assert!(!ms.mipdata_.domains[0].infeasible());
        // lb = min(2, 0) = 0 (no change), ub = max(4, 2) = 4 < 10
        assert_eq!(num_reductions, 1);
        assert_eq!(ms.mipdata_.domains[0].col_upper_[1], 4.0);

        let imp = &*ms.mipdata_.implications;
        assert!(imp.implications[0].computed && imp.implications[1].computed);
        assert_eq!(
            x1_dfp_c_key(&imp.implications[0].implics),
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kLower, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 4.0)
            ]
        );
        assert_eq!(
            x1_dfp_c_key(&imp.implications[1].implics),
            vec![x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0)]
        );
        assert_eq!(imp.num_implications, 3);
        assert!(imp.implications_down.is_empty());
        assert!(imp.implications_up.is_empty());
        assert!(imp.dual_fix_probing_bin_inds_.is_empty());
        assert!(imp.dual_fix_probing_bin_flags_.is_empty());
    }

    #[test]
    fn x1_dfp_c_merge_uses_up_down_on_inactive() {
        // Option on but dual-fix probing inactive (separation and sub-MIP probing):
        // 1.16.0 computeImplications still fills implicationsDown/Up with every implication
        // safe (lines 89-92, 101-104), and runProbing merges those vectors (lines 420-487),
        // never the implication tree. Down (x0 = 0) gives x1 in [2, 4], up gives x1 <= 2:
        // lb = min(2, 0) = 0, ub = max(4, 2) = 4, one reduction.
        // The 1.16 profile has both options on.
        let mut ms = x1_dfp_c_probe_model_keep(true, true);
        assert!(!ms.mipdata_.domains[0].get_dual_fix_probing_active());
        let mut num_reductions: HighsInt = 0;
        let probed = detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_reductions
        ));
        assert!(probed);
        assert!(!ms.mipdata_.domains[0].infeasible());
        assert_eq!(num_reductions, 1);
        assert_eq!(ms.mipdata_.domains[0].col_lower_[1], 0.0);
        assert_eq!(ms.mipdata_.domains[0].col_upper_[1], 4.0);

        let imp = &*ms.mipdata_.implications;
        let mut down = x1_dfp_c_key(&imp.implications_down);
        down.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.total_cmp(&b.2)));
        assert_eq!(
            down,
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kLower, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 4.0)
            ]
        );
        assert_eq!(
            x1_dfp_c_key(&imp.implications_up),
            vec![x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0)]
        );
        assert!(imp.has_probed[0] && imp.has_probed[1]);
        assert!(imp.dual_fix_probing_bin_inds_.is_empty());
    }

    #[test]
    fn x1_dfp_c_merge_uses_up_down_on() {
        // Option on and the global domain's dual-fix probing active: computeImplications
        // fills implicationsDown/Up (binaries erased, sorted by column), runProbing merges
        // those (latest lines 381-420) and clears the tentative state.
        // With the DFP propagation of HiGHS backport-DFP-B, each probe also dual-fixes x1 (cost 1)
        // once its down-locking row r0 is redundant: down (x0 = 0) gives x1 in [2, 2],
        // up (x0 = 1) gives x1 in [0, 0]. The merge then takes the substitution branch
        // (lbDown == ubDown && lbUp == ubUp, HighsImplications.cpp latest 458-468):
        // x1 = 2 + (0 - 2) * x0, one reduction, and the global x1 bounds are untouched.
        let mut ms = x1_dfp_c_probe_model(true);
        assert_eq!(
            ms.mipdata_.implications.dual_fix_probing_bin_flags_.len(),
            3
        );
        {
            let ms_ref = &mut ms;
            detach!(ms_ref.mipdata_.domains[0], |gd| gd
                .dfp_recompute_locks(ms_ref));
        }
        ms.mipdata_.domains[0].set_dual_fix_probing_active(true);
        let mut num_reductions: HighsInt = 0;
        let probed = detach!(ms.mipdata_.implications, |imp| imp.run_probing(
            &mut ms,
            0,
            &mut num_reductions
        ));
        assert!(probed);
        assert!(!ms.mipdata_.domains[0].infeasible());
        assert_eq!(num_reductions, 1);
        assert_eq!(ms.mipdata_.domains[0].col_lower_[1], 0.0);
        assert_eq!(ms.mipdata_.domains[0].col_upper_[1], 10.0);

        let imp = &*ms.mipdata_.implications;
        assert_eq!(imp.substitutions.len(), 1);
        assert_eq!(imp.substitutions[0].substcol, 1);
        assert_eq!(imp.substitutions[0].staycol, 0);
        assert_eq!(imp.substitutions[0].offset, 2.0);
        assert_eq!(imp.substitutions[0].scale, -2.0);
        assert!(imp.colsubstituted[1] != 0);

        let mut down = x1_dfp_c_key(&imp.implications_down);
        down.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.total_cmp(&b.2)));
        assert_eq!(
            down,
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kLower, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 4.0)
            ]
        );
        let mut up = x1_dfp_c_key(&imp.implications_up);
        up.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.total_cmp(&b.2)));
        assert_eq!(
            up,
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 0.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0)
            ]
        );
        // the safe prefix is kept in the 1.15.1 storage, sorted by operator<
        assert_eq!(
            x1_dfp_c_key(&imp.implications[0].implics),
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kLower, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 4.0)
            ]
        );
        assert_eq!(
            x1_dfp_c_key(&imp.implications[1].implics),
            vec![
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 0.0),
                x1_dfp_c_chg(1, HighsBoundType::kUpper, 2.0)
            ]
        );
        assert!(imp.implications[0].computed && imp.implications[1].computed);
        assert!(imp.dual_fix_probing_bin_inds_.is_empty());
        assert!(
            imp.dual_fix_probing_bin_flags_
                .iter()
                .all(|f| f.is_undecided())
        );
    }

    #[test]
    fn x1_dfp_c_tentative_fixing_record() {
        // TentativeFixing::record (HighsImplications.h 84-88): a raised lower bound means the
        // column would be fixed to its upper bound; the first record of a probe wins.
        let mut imp = HighsImplications::default();
        imp.dual_fix_probing_bin_flags_ = vec![HighsImplications_TentativeFixing::default(); 2];
        let up = HighsDomainChange {
            boundval: 1.0,
            column: 1,
            boundtype: HighsBoundType::kLower,
        };
        let down = HighsDomainChange {
            boundval: 0.0,
            column: 1,
            boundtype: HighsBoundType::kUpper,
        };
        imp.record_tentative_cliques(1, &up);
        imp.record_tentative_cliques(1, &down);
        imp.record_tentative_cliques(0, &down);
        assert_eq!(imp.dual_fix_probing_bin_inds_, vec![1]);
        assert_eq!(
            imp.dual_fix_probing_bin_flags_[1],
            HighsImplications_TentativeFixing {
                down_probe: HighsImplications_TentativeFixingDirection::FixLower,
                up_probe: HighsImplications_TentativeFixingDirection::FixUpper,
            }
        );
        imp.clear_tentative_cliques();
        assert!(imp.dual_fix_probing_bin_inds_.is_empty());
        assert!(imp.dual_fix_probing_bin_flags_[1].is_undecided());
    }

    #[test]
    fn x1_agg_a_default_origin_minus_one() {
        assert_eq!(HighsImplications_VarBound::default().origin, -1);
        let imp = HighsImplications::default();
        assert!(imp.row_to_var_bounds.is_empty());
    }

    #[test]
    fn x1_agg_a_row_modified_clears_origin() {
        let mut imp = HighsImplications::default();
        imp.vlbs.resize_with(3, HighsHashTree::default);
        imp.vubs.resize_with(3, HighsHashTree::default);
        imp.row_to_var_bounds.resize_with(4, HighsHashTree::default);

        // vlbs[1] key 2, origin 3
        imp.vlbs[1].insert(HighsHashTableEntry::new(
            2,
            HighsImplications_VarBound {
                coef: -1.0,
                constant: 1.0,
                origin: 3,
            },
        ));
        // vubs[1] key 2, origin 3
        imp.vubs[1].insert(HighsHashTableEntry::new(
            2,
            HighsImplications_VarBound {
                coef: 1.0,
                constant: 0.0,
                origin: 3,
            },
        ));
        // vlbs[0] key 2, origin 0
        imp.vlbs[0].insert(HighsHashTableEntry::new(
            2,
            HighsImplications_VarBound {
                coef: -1.0,
                constant: 1.0,
                origin: 0,
            },
        ));

        // register in row_to_var_bounds[3] the pair (1, 2)
        imp.row_to_var_bounds[3].insert(HighsHashTableEntry::new(1, 2));
        // register in row_to_var_bounds[0] the pair (0, 2)
        imp.row_to_var_bounds[0].insert(HighsHashTableEntry::new(0, 2));

        imp.row_modified(3);

        // vlbs[1] key 2 origin should be -1
        assert_eq!(imp.vlbs[1].find(&2).unwrap().origin, -1);
        // vubs[1] key 2 origin should be -1
        assert_eq!(imp.vubs[1].find(&2).unwrap().origin, -1);
        // row_to_var_bounds[3] should be empty
        let mut count = 0;
        imp.row_to_var_bounds[3].for_each(&mut |_k, _v| {
            count += 1;
        });
        assert_eq!(count, 0);
        // vlbs[0] key 2 origin remains 0 (not touched because origin != row)
        assert_eq!(imp.vlbs[0].find(&2).unwrap().origin, 0);
        // row_to_var_bounds[0] still contains (0, 2)
        assert!(imp.row_to_var_bounds[0].contains(&0));
    }

    #[test]
    fn x1_agg_a_tighten_vlb_inf_safe() {
        let lp = HighsLp::default();
        let mut cb = HighsCallback::default();
        let options = HighsOptions::default();
        let solution = HighsSolution::default();
        let ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        let imp = HighsImplications::default();

        // Finite case
        let mut vlb = HighsImplications_VarBound {
            coef: -2.0,
            constant: 3.0,
            origin: -1,
        };
        let result = imp.tighten_vlb(&ms, &mut vlb, 2.0);
        assert!(result);
        assert_eq!(vlb.coef, -1.0);
        assert_eq!(vlb.constant, 3.0);

        // Inf case: constant = +inf
        let vlb_inf = HighsImplications_VarBound {
            coef: -1.0,
            constant: f64::INFINITY,
            origin: -1,
        };
        // inf-safe minValue is +inf (constant + min(coef,0) = inf + (-1) = inf)
        assert_eq!(var_bound_min_value_inf_safe(&vlb_inf), f64::INFINITY);
        // 1.15.1 min_value would give NaN (HighsCDouble(inf) + (-1) is NaN)
        assert!(var_bound_min_value(&vlb_inf).is_nan());

        // tighten_vlb should return false because redundant_vlb: maxValue_inf_safe = inf + max(-1,0) = inf + 0 = inf
        // inf <= lb + feastol is false for finite lb, so not redundant by that check.
        // But min_value_inf_safe = inf >= lb - epsilon is true for finite lb, so return false.
        let mut vlb_inf_mut = vlb_inf;
        let result = imp.tighten_vlb(&ms, &mut vlb_inf_mut, 2.0);
        assert!(!result);
    }

    /// HiGHS #3264 (mip_aggregate_implications): a VUB stored with an infinite coefficient
    /// (probing implication while the column upper bound was infinite) is tightened by
    /// 1.16.0 `cleanupVub` once the bound is finite (coef 4 for ub 4); the 1.15.1
    /// `HighsCDouble` maxValue is NaN, so the coefficient stays infinite.
    #[test]
    fn x1_agg_cleanup_vub_tightens_inf_coef() {
        for (flag, want) in [(true, 4.0), (false, f64::INFINITY)] {
            let lp = HighsLp::default();
            let mut cb = HighsCallback::default();
            let mut options = HighsOptions::default();
            options.mip_aggregate_implications = flag;
            let solution = HighsSolution::default();
            let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
            let mut imp = HighsImplications::default();
            let mut vub = HighsImplications_VarBound {
                coef: f64::INFINITY,
                constant: 0.0,
                origin: -1,
            };
            let (mut redundant, mut infeasible) = (false, false);
            imp.cleanup_vub(
                &mut ms,
                0,
                1,
                &mut vub,
                4.0,
                &mut redundant,
                &mut infeasible,
                false,
            );
            assert!(!redundant && !infeasible);
            assert_eq!(vub.coef, want, "flag {flag}");
            assert_eq!(vub.constant, 0.0);
        }
    }

    /// HiGHS #3250 (mip_presolve_row_append): HighsImplications.cpp::reindexRows (1.16.0
    /// lines 635-655). Origins move to the new row index, origins of removed rows become
    /// -1, and rowToVarBounds is rebuilt over the new rows only from the surviving origins.
    #[test]
    fn x1b_psr_reindex_rows() {
        let mut imp = HighsImplications::default();
        imp.vlbs.resize_with(3, HighsHashTree::default);
        imp.vubs.resize_with(3, HighsHashTree::default);
        let vb = |origin| HighsImplications_VarBound {
            coef: 1.0,
            constant: 0.0,
            origin,
        };
        // col 1: VLB on col 0 from row 2, VUB on col 2 from row 1 (a cut, removed)
        imp.vlbs[1].insert(HighsHashTableEntry::new(0, vb(2)));
        imp.vubs[1].insert(HighsHashTableEntry::new(2, vb(1)));
        // col 0: VUB on col 2 from row 3, VLB on col 1 without a row origin
        imp.vubs[0].insert(HighsHashTableEntry::new(2, vb(3)));
        imp.vlbs[0].insert(HighsHashTableEntry::new(1, vb(-1)));
        // stale entries of the old row numbering must be dropped
        imp.row_to_var_bounds.resize_with(4, HighsHashTree::default);
        imp.row_to_var_bounds[1].insert(HighsHashTableEntry::new(1, 2));

        // rows 0, 2, 3 survive as 0, 1, 2; row 1 is moved to the cut pool
        imp.reindex_rows(3, &[0, -1, 1, 2]);

        assert_eq!(imp.vlbs[1].find(&0).unwrap().origin, 1);
        assert_eq!(imp.vubs[1].find(&2).unwrap().origin, -1);
        assert_eq!(imp.vubs[0].find(&2).unwrap().origin, 2);
        assert_eq!(imp.vlbs[0].find(&1).unwrap().origin, -1);

        assert_eq!(imp.row_to_var_bounds.len(), 3);
        assert!(imp.row_to_var_bounds[0].empty());
        assert_eq!(imp.row_to_var_bounds[1].find(&1), Some(&0));
        assert_eq!(imp.row_to_var_bounds[2].find(&0), Some(&2));
        let mut n = 0;
        for tree in &imp.row_to_var_bounds {
            tree.for_each(&mut |_, _| n += 1);
        }
        assert_eq!(n, 2);
    }
}
