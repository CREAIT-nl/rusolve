//! HiGHS source: mip/HighsDomain.h (class definitions and header-inline methods)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! Ownership (re-entrant objects and HighsDomain propagators):
//! * A domain owns its propagators (`cutpoolpropagation`, `conflict_pool_propagation`,
//!   `obj_prop_`) as in C++; each propagator names its pool by index into
//!   `ms.mipdata_.cutpools` / `conflictpools` and is registered there under its `DomainKey`.
//! * Pool notifications are replayed lazily by [`HighsDomain::replay_pool_events`], called first
//!   by every `&mut self` method that changes bounds or reads propagator state.
//! * C++ `this == &mipdata_->getDomain()` is [`HighsDomain::is_mip_global_domain`].
//! * Copies register fresh keys (`copy_from` / `op_assign`); `release` is the destructor. A
//!   plain `clone()` duplicates the keys and is only a building block of `copy_from`.

use crate::prelude::*;

mod propagation;
pub use propagation::*;
mod activity;
pub use activity::*;
mod stack;
pub use stack::*;
mod conflict_set;
pub use conflict_set::*;
mod dual_fix_probing;
pub use dual_fix_probing::*;

/// HiGHS: HighsDomain.h::HighsDomain (lines 27-672)
///
/// Removed C++ members:
/// * C++ `HighsMipSolver* mipsolver`: removed (M1); every method that reads it takes `ms`.
///
/// `Clone` copies the propagator keys verbatim (two domains would then share an observer
/// identity); use [`HighsDomain::copy_from`] / [`HighsDomain::op_assign`] for the C++ copy
/// constructor / assignment.
#[derive(Clone, Debug)]
pub struct HighsDomain {
    /// C++: `std::vector<uint8_t> changedcolsflags_`
    pub changedcolsflags_: Vec<u8>,
    /// C++: `std::vector<HighsInt> changedcols_`
    pub changedcols_: Vec<HighsInt>,
    /// C++: `std::vector<std::pair<HighsInt, HighsInt>> propRowNumChangedBounds_`
    pub prop_row_num_changed_bounds_: Vec<(HighsInt, HighsInt)>,
    /// C++: `std::vector<HighsDomainChange> domchgstack_`
    pub domchgstack_: Vec<HighsDomainChange>,
    /// C++: `std::vector<Reason> domchgreason_`
    pub domchgreason_: Vec<HighsDomain_Reason>,
    /// C++: `std::vector<std::pair<double, HighsInt>> prevboundval_`
    pub prevboundval_: Vec<(f64, HighsInt)>,
    /// C++: `std::vector<HighsCDouble> activitymin_`
    pub activitymin_: Vec<HighsCDouble>,
    /// C++: `std::vector<HighsCDouble> activitymax_`
    pub activitymax_: Vec<HighsCDouble>,
    /// C++: `std::vector<HighsInt> activitymininf_`
    pub activitymininf_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> activitymaxinf_`
    pub activitymaxinf_: Vec<HighsInt>,
    /// C++: `std::vector<double> capacityThreshold_`
    pub capacity_threshold_: Vec<f64>,
    /// C++: `std::vector<uint8_t> propagateflags_`
    pub propagateflags_: Vec<u8>,
    /// C++: `std::vector<HighsInt> propagateinds_`
    pub propagateinds_: Vec<HighsInt>,
    /// C++: `ObjectivePropagation objProp_`
    pub obj_prop_: HighsDomain_ObjectivePropagation,
    /// C++: `std::deque<CutpoolPropagation> cutpoolpropagation` (a `Vec`: elements are reached
    /// by index, never by a held reference, so deque address stability is not needed)
    pub cutpoolpropagation: Vec<HighsDomain_CutpoolPropagation>,
    /// C++: `std::deque<ConflictPoolPropagation> conflictPoolPropagation`
    pub conflict_pool_propagation: Vec<HighsDomain_ConflictPoolPropagation>,
    /// C++: `bool infeasible_ = false`
    pub infeasible_: bool,
    /// C++: `Reason infeasible_reason`
    pub infeasible_reason: HighsDomain_Reason,
    /// C++: `HighsInt infeasible_pos`
    pub infeasible_pos: HighsInt,
    /// C++: `std::vector<HighsInt> colLowerPos_`
    pub col_lower_pos_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colUpperPos_`
    pub col_upper_pos_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> branchPos_`
    pub branch_pos_: Vec<HighsInt>,
    /// C++: `HighsHashTable<HighsInt> redundantRows_`
    pub redundant_rows_: HighsHashTable<HighsInt, ()>,
    /// C++: `bool recordRedundantRows_ = false`
    pub record_redundant_rows_: bool,
    /// C++: `std::vector<double> col_lower_`
    pub col_lower_: Vec<f64>,
    /// C++: `std::vector<double> col_upper_`
    pub col_upper_: Vec<f64>,
    /// C++: `bool x1_cdouble_act_contrib` (HiGHS #3184, mip_cdouble_act_contrib)
    pub x1_cdouble_act_contrib: bool,
    /// C++: `bool x1_reduced_stack_fast` (HiGHS #3336, mip_reduced_stack_fast)
    pub x1_reduced_stack_fast: bool,
    /// C++: `DualFixProbingPropagation dualFixProbingPropagation` (HiGHS #3244, mip_dual_fix_probing)
    pub dual_fix_probing_: HighsDomain_DualFixProbingPropagation,
    /// C++: `bool dualFixProbingActive_ = false` (HiGHS #3244, mip_dual_fix_probing)
    pub dual_fix_probing_active: bool,
    /// rusolve (shiftprop): minimum relative improvement of a continuous column's range for
    /// propagation to accept a tightened bound. HiGHS hard-codes 0.3 (HighsDomain.cpp
    /// adjustedLb/adjustedUb); SCIP's default `numerics/boundstreps` is 0.05.
    pub cont_min_improve: f64,
}

impl Default for HighsDomain {
    fn default() -> Self {
        Self {
            changedcolsflags_: Default::default(),
            changedcols_: Default::default(),
            prop_row_num_changed_bounds_: Default::default(),
            domchgstack_: Default::default(),
            domchgreason_: Default::default(),
            prevboundval_: Default::default(),
            activitymin_: Default::default(),
            activitymax_: Default::default(),
            activitymininf_: Default::default(),
            activitymaxinf_: Default::default(),
            capacity_threshold_: Default::default(),
            propagateflags_: Default::default(),
            propagateinds_: Default::default(),
            obj_prop_: Default::default(),
            cutpoolpropagation: Default::default(),
            conflict_pool_propagation: Default::default(),
            infeasible_: false,
            infeasible_reason: Default::default(),
            infeasible_pos: Default::default(),
            col_lower_pos_: Default::default(),
            col_upper_pos_: Default::default(),
            branch_pos_: Default::default(),
            redundant_rows_: Default::default(),
            record_redundant_rows_: false,
            col_lower_: Default::default(),
            col_upper_: Default::default(),
            x1_cdouble_act_contrib: false,
            x1_reduced_stack_fast: false,
            dual_fix_probing_: Default::default(),
            dual_fix_probing_active: false,
            cont_min_improve: 0.3,
        }
    }
}

/// HiGHS: HighsDomain.h::HighsDomain::Reason (lines 29-60)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsDomain_Reason {
    /// C++: `HighsInt type`
    pub r#type: HighsInt,
    /// C++: `HighsInt index`
    pub index: HighsInt,
}
impl HighsDomain_Reason {
    pub const kBranching: HighsInt = -1;
    pub const kUnknown: HighsInt = -2;
    pub const kModelRowUpper: HighsInt = -3;
    pub const kModelRowLower: HighsInt = -4;
    pub const kCliqueTable: HighsInt = -5;
    pub const kConflictingBounds: HighsInt = -6;
    pub const kObjective: HighsInt = -7;
}

/// HiGHS: HighsDomain.h::HighsDomain::ConflictSet (lines 62-157)
///
/// Removed C++ members:
/// * C++: `HighsDomain & localdom`: removed; methods that use it take `localdom: &mut HighsDomain`.
/// * C++: `const HighsDomain & globaldom`: removed; methods that use it take `globaldom: Option<&HighsDomain>`.
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_ConflictSet {
    /// C++: `std::set<LocalDomChg> reasonSideFrontier`
    pub reason_side_frontier: BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
    /// C++: `std::set<LocalDomChg> reconvergenceFrontier`
    pub reconvergence_frontier: BTreeSet<HighsDomain_ConflictSet_LocalDomChg>,
    /// C++: `std::vector<std::set<LocalDomChg>::iterator> resolveQueue`
    pub resolve_queue: Vec<HighsDomain_ConflictSet_LocalDomChg>,
    /// C++: `std::vector<LocalDomChg> resolvedDomainChanges`
    pub resolved_domain_changes: Vec<HighsDomain_ConflictSet_LocalDomChg>,
    /// C++: `std::vector<ResolveCandidate> resolveBuffer`
    pub resolve_buffer: Vec<HighsDomain_ConflictSet_ResolveCandidate>,
}

/// HiGHS: HighsDomain.h::HighsDomain::ConflictSet::LocalDomChg (lines 68-73)
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_ConflictSet_LocalDomChg {
    /// C++: `HighsInt pos`
    pub pos: HighsInt,
    /// C++: `HighsDomainChange domchg`
    pub domchg: HighsDomainChange,
}

/// `std::set`/`std::sort` ordering = C++ `operator<` (equivalence = neither is less).
impl PartialEq for HighsDomain_ConflictSet_LocalDomChg {
    fn eq(&self, other: &Self) -> bool {
        !self.op_lt(other) && !other.op_lt(self)
    }
}
impl Eq for HighsDomain_ConflictSet_LocalDomChg {}
impl PartialOrd for HighsDomain_ConflictSet_LocalDomChg {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HighsDomain_ConflictSet_LocalDomChg {
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

/// HiGHS: HighsDomain.h::HighsDomain::ConflictSet::ResolveCandidate (lines 90-103)
#[derive(Clone, Copy, Default, Debug)]
pub struct HighsDomain_ConflictSet_ResolveCandidate {
    /// C++: `double delta`
    pub delta: f64,
    /// C++: `double baseBound`
    pub base_bound: f64,
    /// C++: `double prio`
    pub prio: f64,
    /// C++: `HighsInt boundPos`
    pub bound_pos: HighsInt,
    /// C++: `HighsInt valuePos`
    pub value_pos: HighsInt,
}

/// `std::set`/`std::sort` ordering = C++ `operator<` (equivalence = neither is less).
impl PartialEq for HighsDomain_ConflictSet_ResolveCandidate {
    fn eq(&self, other: &Self) -> bool {
        !self.op_lt(other) && !other.op_lt(self)
    }
}
impl Eq for HighsDomain_ConflictSet_ResolveCandidate {}
impl PartialOrd for HighsDomain_ConflictSet_ResolveCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HighsDomain_ConflictSet_ResolveCandidate {
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

/// HiGHS: HighsDomain.h::HighsDomain::CutpoolPropagation (lines 159-193)
///
/// Removed C++ members:
/// * C++ `HighsDomain* domain`: removed; methods take the owning domain as a parameter
///   (propagation.rs module note).
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_CutpoolPropagation {
    /// C++: `HighsInt cutpoolindex` (position in the owner's `cutpoolpropagation`; the
    /// `Reason::cut` type, not the pool index)
    pub cutpoolindex: HighsInt,
    /// C++ `HighsCutPool* cutpool`: index of the pool in `ms.mipdata_.cutpools` (= its `index_`).
    pub cutpool: usize,
    /// M3 event-log cursor: absolute index of the next pool event this propagator has not
    /// replayed. `u64::MAX` = not registered (copy made under the parallel lock, C++ skips
    /// `addPropagationDomain`): never replays.
    pub cursor: u64,
    /// M3 observer identity (`PoolEventLog` key) of this propagation object. Not a C++ member:
    /// C++ identifies the observer by its `this` pointer in `propagationDomains`.
    /// PORT NOTE: field added so `HighsCutPool::add_propagation_domain` can register it (cut_pool.rs).
    pub key: crate::mip::DomainKey,
    /// C++: `std::vector<HighsCDouble> activitycuts_`
    pub activitycuts_: Vec<HighsCDouble>,
    /// C++: `std::vector<HighsInt> activitycutsinf_`
    pub activitycutsinf_: Vec<HighsInt>,
    /// C++: `std::vector<uint8_t> propagatecutflags_`
    pub propagatecutflags_: Vec<u8>,
    /// C++: `std::vector<HighsInt> propagatecutinds_`
    pub propagatecutinds_: Vec<HighsInt>,
    /// C++: `std::vector<double> capacityThreshold_`
    pub capacity_threshold_: Vec<f64>,
}

/// HiGHS: HighsDomain.h::HighsDomain::ConflictPoolPropagation (lines 195-236)
///
/// Removed C++ members:
/// * C++ `HighsDomain* domain`: removed; methods take the owning domain as a parameter, and
///   `propagateConflict` is `HighsDomain::propagate_conflict` (propagation.rs module note).
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_ConflictPoolPropagation {
    /// C++: `HighsInt conflictpoolindex` (position in the owner's `conflictPoolPropagation`)
    pub conflictpoolindex: HighsInt,
    /// C++ `HighsConflictPool* conflictpool_`: index of the pool in `ms.mipdata_.conflictpools`.
    pub conflictpool: usize,
    /// M3 event-log cursor (as `HighsDomain_CutpoolPropagation::cursor`, `u64::MAX` = not registered).
    pub cursor: u64,
    /// M3 observer identity, the key `HighsConflictPool::add_propagation_domain` registers.
    /// PORT NOTE: field added because the C++ pointer registration needs an identity the pool's
    /// `PoolEventLog` can key on (`DomainKey`); not in the C++ class.
    pub key: DomainKey,
    /// C++: `std::vector<HighsInt> colLowerWatched_`
    pub col_lower_watched_: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colUpperWatched_`
    pub col_upper_watched_: Vec<HighsInt>,
    /// C++: `std::vector<uint8_t> conflictFlag_`
    pub conflict_flag_: Vec<u8>,
    /// C++: `std::vector<HighsInt> propagateConflictInds_`
    pub propagate_conflict_inds_: Vec<HighsInt>,
    /// C++: `std::vector<WatchedLiteral> watchedLiterals_`
    pub watched_literals_: Vec<HighsDomain_ConflictPoolPropagation_WatchedLiteral>,
}

/// HiGHS: HighsDomain.h::HighsDomain::ConflictPoolPropagation::WatchedLiteral (lines 204-208)
#[derive(Clone, Debug)]
pub struct HighsDomain_ConflictPoolPropagation_WatchedLiteral {
    /// C++: `HighsDomainChange domchg = { 0.0 , - 1 , HighsBoundType :: kLower }`
    pub domchg: HighsDomainChange,
    /// C++: `HighsInt prev = - 1`
    pub prev: HighsInt,
    /// C++: `HighsInt next = - 1`
    pub next: HighsInt,
}

impl Default for HighsDomain_ConflictPoolPropagation_WatchedLiteral {
    fn default() -> Self {
        Self {
            domchg: HighsDomainChange {
                boundval: 0.0,
                column: -1,
                boundtype: HighsBoundType::kLower,
            },
            prev: -1,
            next: -1,
        }
    }
}

/// HiGHS: HighsDomain.h::HighsDomain::ObjectivePropagation (lines 239-298)
///
/// Removed C++ members (the member functions are `HighsDomain::obj_prop_*`, propagation.rs):
/// * C++: `HighsDomain * domain = nullptr`: the owning domain is `self`; `domain != nullptr` is `active`.
/// * C++: `const HighsObjectiveFunction * objFunc`: `ms.mipdata_.objective_function`.
/// * C++: `const double * cost`: `mip_model!(ms).col_cost_`.
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_ObjectivePropagation {
    /// C++ `domain != nullptr` (set by `HighsDomain::obj_prop_init`)
    pub active: bool,
    /// C++: `HighsCDouble objectiveLower`
    pub objective_lower: HighsCDouble,
    /// C++: `HighsInt numInfObjLower`
    pub num_inf_obj_lower: HighsInt,
    /// C++: `double capacityThreshold`
    pub capacity_threshold: f64,
    /// C++: `bool isPropagated`
    pub is_propagated: bool,
    /// C++: `std::vector<ObjectiveContribution> objectiveLowerContributions`
    pub objective_lower_contributions: Vec<HighsDomain_ObjectivePropagation_ObjectiveContribution>,
    /// C++: `std::vector<std::pair<HighsInt, HighsInt>> contributionPartitionSets`
    pub contribution_partition_sets: Vec<(HighsInt, HighsInt)>,
    /// C++: `std::vector<double> propagationConsBuffer`
    pub propagation_cons_buffer: Vec<f64>,
    /// C++: `std::vector<PartitionCliqueData> partitionCliqueData`
    pub partition_clique_data: Vec<HighsDomain_ObjectivePropagation_PartitionCliqueData>,
}

/// HiGHS: HighsDomain.h::HighsDomain::ObjectivePropagation::ObjectiveContribution (lines 248-253)
#[derive(Clone, Default, Debug)]
pub struct HighsDomain_ObjectivePropagation_ObjectiveContribution {
    /// C++: `double contribution`
    pub contribution: f64,
    /// C++: `HighsInt col`
    pub col: HighsInt,
    /// C++: `HighsInt partition`
    pub partition: HighsInt,
    /// C++: `highs::RbTreeLinks<HighsInt> links`
    pub links: RbTreeLinks<HighsInt>,
}

/// HiGHS: HighsDomain.h::HighsDomain::ObjectivePropagation::PartitionCliqueData (lines 260-264)
#[derive(Clone, Copy, Default, Debug, PartialEq)]
pub struct HighsDomain_ObjectivePropagation_PartitionCliqueData {
    /// C++: `double multiplier`
    pub multiplier: f64,
    /// C++: `HighsInt rhs`
    pub rhs: HighsInt,
    /// C++: `bool changed`
    pub changed: bool,
}

impl HighsDomain_Reason {
    /// HiGHS: HighsDomain.h::HighsDomain::Reason::branching (lines 42-42)
    /// C++: `static Reason branching()`
    pub fn branching() -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kBranching,
            index: 0,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::unspecified (lines 43-43)
    /// C++: `static Reason unspecified()`
    pub fn unspecified() -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kUnknown,
            index: 0,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::cliqueTable (lines 44-46)
    /// C++: `static Reason cliqueTable(HighsInt col, HighsInt val)`
    pub fn clique_table(col: HighsInt, val: HighsInt) -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kCliqueTable,
            index: 2 * col + val,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::modelRowUpper (lines 47-49)
    /// C++: `static Reason modelRowUpper(HighsInt row)`
    pub fn model_row_upper(row: HighsInt) -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kModelRowUpper,
            index: row,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::modelRowLower (lines 50-52)
    /// C++: `static Reason modelRowLower(HighsInt row)`
    pub fn model_row_lower(row: HighsInt) -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kModelRowLower,
            index: row,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::cut (lines 53-55)
    /// C++: `static Reason cut(HighsInt cutpool, HighsInt cut)`
    pub fn cut(cutpool: HighsInt, cut: HighsInt) -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: cutpool,
            index: cut,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::conflictingBounds (lines 56-58)
    /// C++: `static Reason conflictingBounds(HighsInt pos)`
    pub fn conflicting_bounds(pos: HighsInt) -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kConflictingBounds,
            index: pos,
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::Reason::objective (lines 59-59)
    /// C++: `static Reason objective()`
    pub fn objective() -> HighsDomain_Reason {
        HighsDomain_Reason {
            r#type: Self::kObjective,
            index: 0,
        }
    }
}

impl HighsDomain_ConflictSet_LocalDomChg {
    /// HiGHS: HighsDomain.h::HighsDomain::ConflictSet::LocalDomChg::operator< (lines 72-72)
    /// C++: `bool operator<(const LocalDomChg & other) const`
    pub fn op_lt(&self, other: &HighsDomain_ConflictSet_LocalDomChg) -> bool {
        self.pos < other.pos
    }
}

impl HighsDomain_ConflictSet_ResolveCandidate {
    /// HiGHS: HighsDomain.h::HighsDomain::ConflictSet::ResolveCandidate::operator< (lines 97-102)
    /// C++: `bool operator<(const ResolveCandidate & other) const`
    pub fn op_lt(&self, other: &HighsDomain_ConflictSet_ResolveCandidate) -> bool {
        if self.prio > other.prio {
            return true;
        }
        if other.prio > self.prio {
            return false;
        }

        self.bound_pos < other.bound_pos
    }
}

impl HighsDomain_ObjectivePropagation {
    /// HiGHS: HighsDomain.h::HighsDomain::ObjectivePropagation::ObjectivePropagation (lines 268-275)
    /// C++: `ObjectivePropagation()`
    pub fn new_2() -> Self {
        // C++ default constructor: objFunc / cost null, objectiveLower 0, numInfObjLower 0,
        // capacityThreshold 0, isPropagated false, domain nullptr (inactive).
        Self::default()
    }

    /// HiGHS: HighsDomain.h::HighsDomain::ObjectivePropagation::isActive (lines 278-278)
    /// C++: `bool isActive() const`
    pub fn is_active(&self) -> bool {
        self.active
    }
}

impl HighsDomain {
    /// HiGHS: HighsDomain.h::HighsDomain::HighsDomain (lines 356-386)
    /// C++: `HighsDomain(const HighsDomain & other)`
    ///
    /// PORT NOTE: signature changed: the propagator copies register with their pools
    /// (`HighsDomain_CutpoolPropagation::copy_from(ms, ..)`), so `ms` is needed. When `other`
    /// lives in `ms` detach it: `detach!(ms.mipdata_.domains[0], |g| HighsDomain::copy_from(ms, g))`.
    /// Members the C++ initializer list omits (`propRowNumChangedBounds_`, `redundantRows_`,
    /// `recordRedundantRows_`) take their defaults.
    pub fn copy_from(ms: &mut HighsMipSolver, other: &HighsDomain) -> Self {
        // C++ copy constructor: `cutpoolpropagation(other.cutpoolpropagation)` copy-constructs
        // every element, in order, and each registers with its pool.
        let mut cutpoolpropagation = Vec::with_capacity(other.cutpoolpropagation.len());
        for prop in other.cutpoolpropagation.iter() {
            cutpoolpropagation.push(HighsDomain_CutpoolPropagation::copy_from(ms, prop));
        }
        let mut conflict_pool_propagation =
            Vec::with_capacity(other.conflict_pool_propagation.len());
        for prop in other.conflict_pool_propagation.iter() {
            conflict_pool_propagation
                .push(HighsDomain_ConflictPoolPropagation::copy_from(ms, prop));
        }
        let mut copy = HighsDomain {
            changedcolsflags_: other.changedcolsflags_.clone(),
            changedcols_: other.changedcols_.clone(),
            prop_row_num_changed_bounds_: Default::default(),
            domchgstack_: other.domchgstack_.clone(),
            domchgreason_: other.domchgreason_.clone(),
            prevboundval_: other.prevboundval_.clone(),
            activitymin_: other.activitymin_.clone(),
            activitymax_: other.activitymax_.clone(),
            activitymininf_: other.activitymininf_.clone(),
            activitymaxinf_: other.activitymaxinf_.clone(),
            capacity_threshold_: other.capacity_threshold_.clone(),
            propagateflags_: other.propagateflags_.clone(),
            propagateinds_: other.propagateinds_.clone(),
            obj_prop_: other.obj_prop_.clone(),
            cutpoolpropagation,
            conflict_pool_propagation,
            infeasible_: other.infeasible_,
            infeasible_reason: other.infeasible_reason,
            infeasible_pos: other.infeasible_pos,
            col_lower_pos_: other.col_lower_pos_.clone(),
            col_upper_pos_: other.col_upper_pos_.clone(),
            branch_pos_: other.branch_pos_.clone(),
            redundant_rows_: Default::default(),
            record_redundant_rows_: false,
            col_lower_: other.col_lower_.clone(),
            col_upper_: other.col_upper_.clone(),
            x1_cdouble_act_contrib: other.x1_cdouble_act_contrib,
            x1_reduced_stack_fast: other.x1_reduced_stack_fast,
            dual_fix_probing_: Default::default(),
            dual_fix_probing_active: false,
            cont_min_improve: other.cont_min_improve,
        };
        if ms.parallel_lock_active() {
            copy.replay_unregistered_backlog(ms, other);
        }
        copy
    }

    /// Parallel-lock copy (parallel glue, no C++ counterpart): bring the pool-0 propagators the
    /// copy did not register (`cursor == u64::MAX`, C++ skips the registration under
    /// `parallelLockActive()`) up to date with the pool's current log.
    ///
    /// C++ delivers every pool notification synchronously, so the source domain's pool-0
    /// activities are always current and the unregistered copy never needs another event
    /// (pool 0 does not change while the tasks run). Here the source replays lazily
    /// (`other`'s cursor may trail the log: pool 0 grew in a serial phase and `other`, a
    /// worker's global domain, has not run a method since), so a copy that inherited the
    /// arrays and then dropped the cursor stayed stale for good and the first bound change
    /// indexed past its arrays. Replaying `other`'s backlog on the copy, without logging an
    /// advance (the copy has no registered key), restores the C++ state.
    fn replay_unregistered_backlog(&mut self, ms: &HighsMipSolver, other: &HighsDomain) {
        if !self.cutpoolpropagation.is_empty() {
            let mut props = std::mem::take(&mut self.cutpoolpropagation);
            for (prop, src) in props.iter_mut().zip(other.cutpoolpropagation.iter()) {
                if prop.cursor != u64::MAX || src.cursor == u64::MAX {
                    continue;
                }
                let log = &ms.mipdata_.cutpools[prop.cutpool].events;
                if src.cursor >= log.end() {
                    continue;
                }
                for event in log.pending(src.cursor) {
                    match event {
                        PoolEvent::CutAdded {
                            cut,
                            propagate,
                            inds,
                            vals,
                            rhs,
                        } => prop.cut_added(ms, self, *cut, *propagate, inds, vals, *rhs),
                        PoolEvent::CutDeleted { cut, only_for_prop } => {
                            prop.cut_deleted(ms, self, *cut, *only_for_prop)
                        }
                        _ => unreachable!("conflict event in a cut pool log"),
                    }
                }
            }
            self.cutpoolpropagation = props;
        }
        if !self.conflict_pool_propagation.is_empty() {
            let mut props = std::mem::take(&mut self.conflict_pool_propagation);
            for (prop, src) in props.iter_mut().zip(other.conflict_pool_propagation.iter()) {
                if prop.cursor != u64::MAX || src.cursor == u64::MAX {
                    continue;
                }
                let log = &ms.mipdata_.conflictpools[prop.conflictpool].events;
                if src.cursor >= log.end() {
                    continue;
                }
                for event in log.pending(src.cursor) {
                    match event {
                        PoolEvent::ConflictAdded { conflict, entries } => prop.conflict_added(
                            self,
                            *conflict,
                            entries,
                            ms.options_mip_.h16_watch_link,
                        ),
                        PoolEvent::ConflictDeleted { conflict } => {
                            prop.conflict_deleted(*conflict, ms.options_mip_.h16_watch_link)
                        }
                        _ => unreachable!("cut event in a conflict pool log"),
                    }
                }
            }
            self.conflict_pool_propagation = props;
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::operator= (lines 388-418)
    /// C++: `HighsDomain & operator=(const HighsDomain & other)`
    ///
    /// PORT NOTE: `cutpoolpropagation = other.cutpoolpropagation` is `std::deque::operator=`:
    /// element-wise `op_assign(ms, ..)` for the common prefix, `copy_from(ms, ..)` for extra
    /// elements, `release(ms)` for surplus ones (in that order). Members the C++ leaves
    /// untouched (`propRowNumChangedBounds_`, `infeasible_pos`, `redundantRows_`,
    /// `recordRedundantRows_`) keep their values. Same detach idiom as `copy_from` when `other`
    /// lives in `ms` (`localdom = mipdata_->getDomain()`).
    pub fn op_assign(&mut self, ms: &mut HighsMipSolver, other: &HighsDomain) {
        self.changedcolsflags_ = other.changedcolsflags_.clone();
        self.changedcols_ = other.changedcols_.clone();
        self.domchgstack_ = other.domchgstack_.clone();
        self.domchgreason_ = other.domchgreason_.clone();
        self.prevboundval_ = other.prevboundval_.clone();
        self.activitymin_ = other.activitymin_.clone();
        self.activitymax_ = other.activitymax_.clone();
        self.activitymininf_ = other.activitymininf_.clone();
        self.activitymaxinf_ = other.activitymaxinf_.clone();
        self.capacity_threshold_ = other.capacity_threshold_.clone();
        self.propagateflags_ = other.propagateflags_.clone();
        self.propagateinds_ = other.propagateinds_.clone();
        self.obj_prop_ = other.obj_prop_.clone();
        // C++ `cutpoolpropagation = other.cutpoolpropagation` is `std::deque::operator=`:
        // element-wise assignment of the common prefix, then copy-construction of the extra
        // elements, or destruction (release) of the surplus ones.
        let n_self = self.cutpoolpropagation.len();
        let n_other = other.cutpoolpropagation.len();
        let common = std::cmp::min(n_self, n_other);
        for i in 0..common {
            self.cutpoolpropagation[i].op_assign(ms, &other.cutpoolpropagation[i]);
        }
        if n_other > n_self {
            for prop in other.cutpoolpropagation[common..].iter() {
                self.cutpoolpropagation
                    .push(HighsDomain_CutpoolPropagation::copy_from(ms, prop));
            }
        } else {
            let surplus: Vec<HighsDomain_CutpoolPropagation> =
                self.cutpoolpropagation.drain(common..).collect();
            for mut prop in surplus {
                prop.release(ms);
            }
        }
        // Same for `conflictPoolPropagation`.
        let n_self = self.conflict_pool_propagation.len();
        let n_other = other.conflict_pool_propagation.len();
        let common = std::cmp::min(n_self, n_other);
        for i in 0..common {
            self.conflict_pool_propagation[i].op_assign(ms, &other.conflict_pool_propagation[i]);
        }
        if n_other > n_self {
            for prop in other.conflict_pool_propagation[common..].iter() {
                self.conflict_pool_propagation
                    .push(HighsDomain_ConflictPoolPropagation::copy_from(ms, prop));
            }
        } else {
            let surplus: Vec<HighsDomain_ConflictPoolPropagation> =
                self.conflict_pool_propagation.drain(common..).collect();
            for mut prop in surplus {
                prop.release(ms);
            }
        }
        self.infeasible_ = other.infeasible_;
        self.infeasible_reason = other.infeasible_reason;
        // C++ `operator=` does not copy `infeasible_pos` (nor `propRowNumChangedBounds_`,
        // `redundantRows_`, `recordRedundantRows_`): those keep their values.
        self.col_lower_pos_ = other.col_lower_pos_.clone();
        self.col_upper_pos_ = other.col_upper_pos_.clone();
        self.branch_pos_ = other.branch_pos_.clone();
        self.col_lower_ = other.col_lower_.clone();
        self.col_upper_ = other.col_upper_.clone();
        self.x1_cdouble_act_contrib = other.x1_cdouble_act_contrib;
        self.x1_reduced_stack_fast = other.x1_reduced_stack_fast;
        // C++ operator=: dualFixProbingActive_ = false; endProbing() (the DFP struct is not copied)
        self.dual_fix_probing_active = false;
        self.dual_fix_probing_.end_probing();
    }

    /// C++ `~HighsDomain()` (implicit; runs `~CutpoolPropagation` / `~ConflictPoolPropagation`):
    /// unregister every propagator from its pool. Not a C++ member function.
    ///
    /// PORT NOTE: added (M3): call before dropping a domain that registered with pools.
    pub fn release(&mut self, ms: &mut HighsMipSolver) {
        for prop in self.cutpoolpropagation.iter_mut() {
            prop.release(ms);
        }
        for prop in self.conflict_pool_propagation.iter_mut() {
            prop.release(ms);
        }
    }

    /// C++ `this == &mipdata_->getDomain()` (HighsDomain.cpp:450, 490; also HighsSeparation,
    /// HighsLpRelaxation). Not a C++ member function.
    ///
    /// PORT NOTE: added: identity through the slot's stable address, valid while
    /// `domains[0]` is attached, detached (this method running on it) or lent.
    pub fn is_mip_global_domain(&self, ms: &HighsMipSolver) -> bool {
        ms.mipdata_.domains[0].holds(self)
    }

    /// M3: replay, in log order, every pool event this domain's propagators have not seen
    /// (C++ delivers them synchronously: `cutAdded`, `cutDeleted`, `conflictAdded`,
    /// `conflictDeleted`). Not a C++ member function.
    ///
    /// PORT NOTE: added. Called first by every `&mut self` method that changes bounds or
    /// reads propagator state (`change_bound`, `backtrack*`, `set_domain_change_stack*`,
    /// `propagate`, `conflict_analysis*`, `get_min_cut_activity`). Events of different pools
    /// commute (each propagator only reads its own state and the unchanged domain bounds), so
    /// replaying pool by pool equals the C++ interleaving. The pools must be attached.
    /// PORT NOTE (perf): true iff [`HighsDomain::replay_pool_events`] has something to do,
    /// i.e. some pool log has events past this domain's cursor. Lets `&`-callers skip the
    /// detach-and-replay round trip when it would be a no-op.
    pub fn has_pending_pool_events(&self, ms: &HighsMipSolver) -> bool {
        self.cutpoolpropagation
            .iter()
            .any(|prop| prop.cursor < ms.mipdata_.cutpools[prop.cutpool].events.end())
            || self
                .conflict_pool_propagation
                .iter()
                .any(|prop| prop.cursor < ms.mipdata_.conflictpools[prop.conflictpool].events.end())
    }

    pub fn replay_pool_events(&mut self, ms: &HighsMipSolver) {
        // Nothing pending: return without touching `self` (the global domain is shared
        // by concurrent worker tasks, see `HighsMipSolver::pin_shared_slots`).
        if !self.has_pending_pool_events(ms) {
            return;
        }
        if !self.cutpoolpropagation.is_empty() {
            let mut props = std::mem::take(&mut self.cutpoolpropagation);
            for prop in props.iter_mut() {
                let log = &ms.mipdata_.cutpools[prop.cutpool].events;
                let end = log.end();
                if prop.cursor >= end {
                    continue;
                }
                for event in log.pending(prop.cursor) {
                    match event {
                        PoolEvent::CutAdded {
                            cut,
                            propagate,
                            inds,
                            vals,
                            rhs,
                        } => prop.cut_added(ms, self, *cut, *propagate, inds, vals, *rhs),
                        PoolEvent::CutDeleted { cut, only_for_prop } => {
                            prop.cut_deleted(ms, self, *cut, *only_for_prop)
                        }
                        _ => unreachable!("conflict event in a cut pool log"),
                    }
                }
                prop.cursor = end;
                log.advance(prop.key, end);
            }
            self.cutpoolpropagation = props;
        }
        if !self.conflict_pool_propagation.is_empty() {
            let mut props = std::mem::take(&mut self.conflict_pool_propagation);
            for prop in props.iter_mut() {
                let log = &ms.mipdata_.conflictpools[prop.conflictpool].events;
                let end = log.end();
                if prop.cursor >= end {
                    continue;
                }
                for event in log.pending(prop.cursor) {
                    match event {
                        // HiGHS 1.16 #3571 (h16_watch_link): read at each replay, passed through (1.16.0 HighsDomain.cpp)
                        PoolEvent::ConflictAdded { conflict, entries } => prop.conflict_added(
                            self,
                            *conflict,
                            entries,
                            ms.options_mip_.h16_watch_link,
                        ),
                        PoolEvent::ConflictDeleted { conflict } => {
                            prop.conflict_deleted(*conflict, ms.options_mip_.h16_watch_link)
                        }
                        _ => unreachable!("cut event in a conflict pool log"),
                    }
                }
                prop.cursor = end;
                log.advance(prop.key, end);
            }
            self.conflict_pool_propagation = props;
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getChangedCols (lines 442-442)
    /// C++: `const std::vector<HighsInt> & getChangedCols() const`
    pub fn get_changed_cols(&self) -> &[HighsInt] {
        &self.changedcols_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::clearChangedCols (lines 448-451)
    /// C++: `void clearChangedCols()`
    pub fn clear_changed_cols(&mut self) {
        for &i in self.changedcols_.iter() {
            self.changedcolsflags_[i as usize] = 0;
        }
        self.changedcols_.clear();
    }

    /// HiGHS: HighsDomain.h::HighsDomain::removeContinuousChangedCols (lines 453-461)
    /// C++: `void removeContinuousChangedCols()`
    pub fn remove_continuous_changed_cols(&mut self, ms: &HighsMipSolver) {
        for &i in self.changedcols_.iter() {
            self.changedcolsflags_[i as usize] = ms.is_col_integral(i) as u8;
        }

        // `erase(remove_if(p))` with `p = ![&](i) { return !isChangedCol(i); }`: keep the
        // columns whose flag is set. The flags are not written during the pass.
        let flags = &self.changedcolsflags_;
        self.changedcols_.retain(|&i| flags[i as usize] != 0);
    }

    /// HiGHS: HighsDomain.h::HighsDomain::clearChangedCols (lines 463-468)
    /// C++: `void clearChangedCols(size_t start)`
    pub fn clear_changed_cols_2(&mut self, start: usize) {
        for i in start..self.changedcols_.len() {
            self.changedcolsflags_[self.changedcols_[i] as usize] = 0;
        }

        self.changedcols_.truncate(start);
    }

    /// HiGHS: HighsDomain.h::HighsDomain::isChangedCol (lines 470-470)
    /// C++: `bool isChangedCol(HighsInt col) const`
    pub fn is_changed_col(&self, col: HighsInt) -> bool {
        self.changedcolsflags_[col as usize] != 0
    }

    /// HiGHS: HighsDomain.h::HighsDomain::isActive (lines 474-478)
    /// C++: `bool isActive(const HighsDomainChange & domchg) const`
    pub fn is_active(&self, domchg: &HighsDomainChange) -> bool {
        if domchg.boundtype == HighsBoundType::kLower {
            domchg.boundval <= self.col_lower_[domchg.column as usize]
        } else {
            domchg.boundval >= self.col_upper_[domchg.column as usize]
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::setupObjectivePropagation (lines 482-482)
    /// C++: `void setupObjectivePropagation()`
    // PORT NOTE: signature changed: `objProp_ = ObjectivePropagation(this)` is `self.obj_prop_init(ms)`.
    pub fn setup_objective_propagation(&mut self, ms: &HighsMipSolver) {
        self.obj_prop_init(ms);
    }

    /// HiGHS: HighsDomain.h::HighsDomain::markInfeasible (lines 486-490)
    /// C++: `void markInfeasible(Reason reason = Reason :: unspecified ( ))`
    /// C++ default arguments: `reason = Reason :: unspecified ( )`
    pub fn mark_infeasible(&mut self, reason: HighsDomain_Reason) {
        self.infeasible_ = true;
        self.infeasible_pos = self.domchgstack_.len() as HighsInt;
        self.infeasible_reason = reason;
    }

    /// HiGHS: HighsDomain.h::HighsDomain::infeasible (lines 492-492)
    /// C++: `bool infeasible() const`
    pub fn infeasible(&self) -> bool {
        self.infeasible_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::changeBound (lines 497-500)
    /// C++: `void changeBound(HighsBoundType boundtype, HighsInt col, double boundval, Reason reason = Reason :: branching ( ))`
    /// C++ default arguments: `reason = Reason :: branching ( )`
    ///
    /// PORT NOTE: signature changed because the overload forwards to
    /// `changeBound(HighsDomainChange, Reason)`, which takes `ms` (activity.rs).
    pub fn change_bound_2(
        &mut self,
        ms: &mut HighsMipSolver,
        boundtype: HighsBoundType,
        col: HighsInt,
        boundval: f64,
        reason: HighsDomain_Reason,
    ) {
        self.change_bound(
            ms,
            HighsDomainChange {
                boundval,
                column: col,
                boundtype,
            },
            reason,
        );
    }

    /// HiGHS: HighsDomain.h::HighsDomain::fixCol (lines 505-516)
    /// C++: `void fixCol(HighsInt col, double val, Reason reason = Reason :: unspecified ( ))`
    /// C++ default arguments: `reason = Reason :: unspecified ( )`
    pub fn fix_col(
        &mut self,
        ms: &mut HighsMipSolver,
        col: HighsInt,
        val: f64,
        reason: HighsDomain_Reason,
    ) {
        // C++ `if (kAllowDeveloperAssert) assert(infeasible_ == 0);`: the flag is false.
        if kAllowDeveloperAssert {
            debug_assert!(!self.infeasible_);
        }
        if self.col_lower_[col as usize] < val {
            self.change_bound(
                ms,
                HighsDomainChange {
                    boundval: val,
                    column: col,
                    boundtype: HighsBoundType::kLower,
                },
                reason,
            );
            if !self.infeasible_ {
                self.propagate(ms);
            }
        }

        if !self.infeasible_ && self.col_upper_[col as usize] > val {
            self.change_bound(
                ms,
                HighsDomainChange {
                    boundval: val,
                    column: col,
                    boundtype: HighsBoundType::kUpper,
                },
                reason,
            );
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getBranchingPositions (lines 522-524)
    /// C++: `const std::vector<HighsInt> & getBranchingPositions() const`
    pub fn get_branching_positions(&self) -> &[HighsInt] {
        &self.branch_pos_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getPreviousBounds (lines 526-528)
    /// C++: `const std::vector<std::pair<double, HighsInt>> & getPreviousBounds() const`
    pub fn get_previous_bounds(&self) -> &[(f64, HighsInt)] {
        &self.prevboundval_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getDomainChangeStack (lines 530-532)
    /// C++: `const std::vector<HighsDomainChange> & getDomainChangeStack() const`
    pub fn get_domain_change_stack(&self) -> &[HighsDomainChange] {
        &self.domchgstack_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getDomainChangeReason (lines 534-536)
    /// C++: `const std::vector<Reason> & getDomainChangeReason() const`
    pub fn get_domain_change_reason(&self) -> &[HighsDomain_Reason] {
        &self.domchgreason_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getObjectiveLowerBound (lines 538-543)
    /// C++: `double getObjectiveLowerBound() const`
    pub fn get_objective_lower_bound(&self) -> f64 {
        if self.obj_prop_.is_active() && self.obj_prop_.num_inf_obj_lower == 0 {
            return f64::from(self.obj_prop_.objective_lower);
        }

        -kHighsInf
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getCutoffConstraint (lines 545-549)
    /// C++: `void getCutoffConstraint(const double *& vals, const HighsInt *& inds, HighsInt & len, double & rhs)`
    // PORT NOTE: signature changed: forwards to `obj_prop_get_propagation_constraint`, whose
    // slices borrow `ms` (objective data) or `self` (the buffer).
    pub fn get_cutoff_constraint<'s>(
        &'s mut self,
        ms: &'s HighsMipSolver,
        vals: &mut &'s [f64],
        inds: &mut &'s [HighsInt],
        len: &mut HighsInt,
        rhs: &mut f64,
    ) {
        // C++ `objProp_.getPropagationConstraint(domchgstack_.size(), vals, inds, len, rhs)`
        // with the default `domchgCol = -1`.
        let stack_size = self.domchgstack_.len() as HighsInt;
        self.obj_prop_get_propagation_constraint(ms, stack_size, vals, inds, len, rhs, -1);
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getNumDomainChanges (lines 551-551)
    /// C++: `HighsInt getNumDomainChanges() const`
    pub fn get_num_domain_changes(&self) -> HighsInt {
        self.domchgstack_.len() as HighsInt
    }

    /// HiGHS: HighsDomain.h::HighsDomain::colBoundsAreGlobal (lines 553-555)
    /// C++: `bool colBoundsAreGlobal(HighsInt col) const`
    pub fn col_bounds_are_global(&self, col: HighsInt) -> bool {
        self.col_lower_pos_[col as usize] == -1 && self.col_upper_pos_[col as usize] == -1
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getBranchDepth (lines 557-557)
    /// C++: `HighsInt getBranchDepth() const`
    pub fn get_branch_depth(&self) -> HighsInt {
        self.branch_pos_.len() as HighsInt
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getReducedDomainChangeStack (lines 559-590)
    /// HiGHS latest (#3336, a59c1b9f4): HighsDomain.h::HighsDomain::getReducedDomainChangeStack (lines 573-630)
    /// C++: `std::vector<HighsDomainChange> getReducedDomainChangeStack(std::vector<HighsInt> & branchingPositions) const`
    pub fn get_reduced_domain_change_stack(
        &self,
        branching_positions: &mut Vec<HighsInt>,
    ) -> Vec<HighsDomainChange> {
        // HiGHS #3336 (mip_reduced_stack_fast): HighsDomain.h lines 575-630 at a59c1b9f4
        if self.x1_reduced_stack_fast {
            if self.domchgstack_.is_empty() {
                return Vec::new();
            }
            let stack_size = self.domchgstack_.len() as HighsInt;
            let num_cols = self.col_lower_.len() as HighsInt;
            let mut log_num_cols: HighsInt = 0;
            let mut v = num_cols - 1;
            while v > 0 {
                log_num_cols += 1;
                v >>= 1;
            }
            if (stack_size - 1) / num_cols <= log_num_cols {
                // Stack smaller than n * (1 + ceil(log2(n)))
                let mut reducedstack: Vec<HighsDomainChange> =
                    Vec::with_capacity(stack_size as usize);
                branching_positions.reserve(self.branch_pos_.len());
                for i in 0..stack_size {
                    let iu = i as usize;
                    let domchg = &self.domchgstack_[iu];
                    // keep only the tightest bound change for each variable
                    if (domchg.boundtype == HighsBoundType::kLower
                        && self.col_lower_pos_[domchg.column as usize] != i)
                        || (domchg.boundtype == HighsBoundType::kUpper
                            && self.col_upper_pos_[domchg.column as usize] != i)
                    {
                        continue;
                    }
                    if self.domchgreason_[iu].r#type == HighsDomain_Reason::kBranching {
                        branching_positions.push(reducedstack.len() as HighsInt);
                    } else {
                        let mut k = i;
                        while self.prevboundval_[k as usize].1 != -1 {
                            k = self.prevboundval_[k as usize].1;
                            if self.domchgreason_[k as usize].r#type
                                == HighsDomain_Reason::kBranching
                            {
                                branching_positions.push(reducedstack.len() as HighsInt);
                                break;
                            }
                        }
                    }
                    reducedstack.push(self.domchgstack_[iu].clone());
                }
                reducedstack.shrink_to_fit();
                return reducedstack;
            } else {
                // Stack is large. Avoid iterating over entire stack to
                // hedge against worst-case performance
                let mut positions: Vec<HighsInt> = Vec::with_capacity(2 * num_cols as usize);
                for col in 0..num_cols {
                    if self.col_lower_pos_[col as usize] != -1 {
                        positions.push(self.col_lower_pos_[col as usize]);
                    }
                    if self.col_upper_pos_[col as usize] != -1 {
                        positions.push(self.col_upper_pos_[col as usize]);
                    }
                }
                pdqsort(&mut positions);
                let mut reducedstack: Vec<HighsDomainChange> = Vec::with_capacity(positions.len());
                branching_positions.reserve(cpp_min(
                    positions.len() as HighsInt,
                    self.branch_pos_.len() as HighsInt,
                ) as usize);
                for &i in positions.iter() {
                    let iu = i as usize;
                    if self.domchgreason_[iu].r#type == HighsDomain_Reason::kBranching {
                        branching_positions.push(reducedstack.len() as HighsInt);
                    } else {
                        let mut k = i;
                        while self.prevboundval_[k as usize].1 != -1 {
                            k = self.prevboundval_[k as usize].1;
                            if self.domchgreason_[k as usize].r#type
                                == HighsDomain_Reason::kBranching
                            {
                                branching_positions.push(reducedstack.len() as HighsInt);
                                break;
                            }
                        }
                    }
                    reducedstack.push(self.domchgstack_[iu].clone());
                }
                reducedstack.shrink_to_fit();
                return reducedstack;
            }
        }

        // 1.15.1 path (off)
        let mut reducedstack: Vec<HighsDomainChange> = Vec::with_capacity(self.domchgstack_.len());
        for i in 0..self.domchgstack_.len() as HighsInt {
            let iu = i as usize;
            let domchg = &self.domchgstack_[iu];
            // keep only the tightest bound change for each variable
            if (domchg.boundtype == HighsBoundType::kLower
                && self.col_lower_pos_[domchg.column as usize] != i)
                || (domchg.boundtype == HighsBoundType::kUpper
                    && self.col_upper_pos_[domchg.column as usize] != i)
            {
                continue;
            }

            if self.domchgreason_[iu].r#type == HighsDomain_Reason::kBranching {
                branching_positions.push(reducedstack.len() as HighsInt);
            } else {
                let mut k = i;
                while self.prevboundval_[k as usize].1 != -1 {
                    k = self.prevboundval_[k as usize].1;
                    if self.domchgreason_[k as usize].r#type == HighsDomain_Reason::kBranching {
                        branching_positions.push(reducedstack.len() as HighsInt);
                        break;
                    }
                }
            }

            reducedstack.push(self.domchgstack_[iu].clone());
        }

        reducedstack
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getMinActivity (lines 622-624)
    /// C++: `double getMinActivity(HighsInt row) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ `getMinActivity`
    // header inline; `HighsTransformedLp::transform` reads it per slack nonzero.
    #[inline]
    pub fn get_min_activity(&self, row: HighsInt) -> f64 {
        if self.activitymininf_[row as usize] == 0 {
            f64::from(self.activitymin_[row as usize])
        } else {
            -kHighsInf
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getMaxActivity (lines 626-628)
    /// C++: `double getMaxActivity(HighsInt row) const`
    // PORT NOTE (perf): `#[inline]`, as the C++ `getMaxActivity`.
    #[inline]
    pub fn get_max_activity(&self, row: HighsInt) -> f64 {
        if self.activitymaxinf_[row as usize] == 0 {
            f64::from(self.activitymax_[row as usize])
        } else {
            kHighsInf
        }
    }

    /// HiGHS: HighsDomain.h::HighsDomain::isBinary (lines 632-635)
    /// C++: `bool isBinary(HighsInt col) const`
    // PORT NOTE: signature changed because C++ isBinary is const; the solver is only read.
    pub fn is_binary(&self, ms: &HighsMipSolver, col: HighsInt) -> bool {
        ms.is_col_integral(col)
            && self.col_lower_[col as usize] == 0.0
            && self.col_upper_[col as usize] == 1.0
    }

    /// HiGHS: HighsDomain.h::HighsDomain::isGlobalBinary (lines 637-641)
    /// C++: `bool isGlobalBinary(HighsInt col) const`
    pub fn is_global_binary(&self, ms: &HighsMipSolver, col: HighsInt) -> bool {
        ms.is_col_integral(col)
            && mip_model!(ms).col_lower_[col as usize] == 0.0
            && mip_model!(ms).col_upper_[col as usize] == 1.0
    }

    /// HiGHS: HighsDomain.h::HighsDomain::variableType (lines 643-645)
    /// C++: `HighsVarType variableType(HighsInt col) const`
    pub fn variable_type(&self, ms: &HighsMipSolver, col: HighsInt) -> HighsVarType {
        ms.variable_type_2(col)
    }

    /// HiGHS: HighsDomain.h::HighsDomain::isFixed (lines 647-649)
    /// C++: `bool isFixed(HighsInt col) const`
    pub fn is_fixed(&self, col: HighsInt) -> bool {
        self.col_lower_[col as usize] == self.col_upper_[col as usize]
    }

    /// HiGHS latest (#3286, 843e9fdb6): HighsDomain.h::HighsDomain::isFixedToVal (lines 666-668)
    /// C++: `bool isFixedToVal(HighsInt col, HighsInt val) const`
    pub fn is_fixed_to_val(&self, col: HighsInt, val: HighsInt) -> bool {
        self.col_lower_[col as usize] == val as f64 && self.col_upper_[col as usize] == val as f64
    }

    /// HiGHS: HighsDomain.h::HighsDomain::numModelNonzeros (lines 657-657)
    /// C++: `HighsInt numModelNonzeros() const`
    pub fn num_model_nonzeros(&self, ms: &HighsMipSolver) -> HighsInt {
        ms.num_nonzero()
    }

    /// HiGHS: HighsDomain.h::HighsDomain::inSubmip (lines 659-659)
    /// C++: `bool inSubmip() const`
    pub fn in_submip(&self, ms: &HighsMipSolver) -> bool {
        ms.submip
    }

    /// HiGHS: HighsDomain.h::HighsDomain::clearRedundantRows (lines 661-661)
    /// C++: `void clearRedundantRows()`
    pub fn clear_redundant_rows(&mut self) {
        self.redundant_rows_.clear()
    }

    /// HiGHS: HighsDomain.h::HighsDomain::getRedundantRows (lines 663-665)
    /// C++: `const HighsHashTable<HighsInt> & getRedundantRows() const`
    pub fn get_redundant_rows(&self) -> &HighsHashTable<HighsInt, ()> {
        &self.redundant_rows_
    }

    /// HiGHS: HighsDomain.h::HighsDomain::setRecordRedundantRows (lines 669-669)
    /// C++: `void setRecordRedundantRows(bool val)`
    pub fn set_record_redundant_rows(&mut self, val: bool) {
        self.record_redundant_rows_ = val;
    }
}

/// HiGHS: HighsDomain.cpp::HighsDomain::ObjectivePropagation::ObjectiveContributionTree (lines 641-670)
///
/// CRTP `CacheMinRbTree` view (key `std::pair<double, HighsInt>`, link
/// `HighsInt`): C++ holds `std::vector<ObjectiveContribution>& nodes =
/// objProp->objectiveLowerContributions` and the base holds
/// `contributionPartitionSets[partition].first` / `.second` as root / first.
pub struct HighsDomain_ObjectivePropagation_ObjectiveContributionTree<'a> {
    pub nodes: &'a mut Vec<HighsDomain_ObjectivePropagation_ObjectiveContribution>,
    pub root: &'a mut HighsInt,
    pub first: &'a mut HighsInt,
}

impl<'a> HighsDomain_ObjectivePropagation_ObjectiveContributionTree<'a> {
    /// HiGHS: HighsDomain.cpp::ObjectiveContributionTree::ObjectiveContributionTree (lines 654-658)
    pub fn new(obj_prop: &'a mut HighsDomain_ObjectivePropagation, partition: HighsInt) -> Self {
        let (root, first) = &mut obj_prop.contribution_partition_sets[partition as usize];
        Self {
            nodes: &mut obj_prop.objective_lower_contributions,
            root,
            first,
        }
    }
}

impl RbTree for HighsDomain_ObjectivePropagation_ObjectiveContributionTree<'_> {
    type KeyType = (f64, HighsInt);
    type LinkType = HighsInt;
    fn root_node(&self) -> HighsInt {
        *self.root
    }
    fn set_root_node(&mut self, root: HighsInt) {
        *self.root = root;
    }
    fn get_rb_tree_links(&self, node: HighsInt) -> &RbTreeLinks<HighsInt> {
        &self.nodes[node as usize].links
    }
    fn get_rb_tree_links_mut(&mut self, node: HighsInt) -> &mut RbTreeLinks<HighsInt> {
        &mut self.nodes[node as usize].links
    }
    /// HiGHS: HighsDomain.cpp::ObjectiveContributionTree::getKey (lines 668-670)
    fn get_key(&self, node: HighsInt) -> (f64, HighsInt) {
        (
            -self.nodes[node as usize].contribution,
            self.nodes[node as usize].col,
        )
    }
}

impl CacheMinRbTree for HighsDomain_ObjectivePropagation_ObjectiveContributionTree<'_> {
    fn first_cached_ref(&self) -> HighsInt {
        *self.first
    }
    fn set_first_cached(&mut self, first: HighsInt) {
        *self.first = first;
    }
}
