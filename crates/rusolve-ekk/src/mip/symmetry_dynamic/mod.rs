//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! constants, sym_* comparisons, bound-tightening wrappers, path-node and
//! symmetry-handling data structs for the dynamic symmetry propagation behind
//! `mip_symmetry_dynamic`.

use crate::prelude::*;
use std::sync::Arc;

pub mod debugsol;
pub mod hook;
pub mod lexred;
pub mod lexred_data;
pub mod lexred_peek;
pub mod lexred_prop;
pub mod orbitopal;
pub mod orbitopal_faces;
pub mod orbitopal_order;
pub mod orbitopal_path;
pub mod setup;
#[cfg(test)]
mod tests_e2e;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING_NONE: HighsInt = 0;
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING_FIRST: HighsInt = 1;
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING_LAST: HighsInt = 2;
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING_CENTRE: HighsInt = 3;
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING_MEDIAN: HighsInt = 4;

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_ColumnOrdering (lines 48-56)
pub(crate) const S3_COLUMN_ORDERING: HighsInt = S3_COLUMN_ORDERING_MEDIAN;

/// Maximum number of propagation rounds.
pub(crate) const S3_MAX_PROP_ROUNDS: HighsInt = 100;

// ---------------------------------------------------------------------------
// sym_* comparisons
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry.c::SCIPsymEQ (lines 1529-1559)
pub(crate) fn sym_eq(a: f64, b: f64, eps: f64, inf: f64) -> bool {
    let inf1 = a >= inf;
    let inf2 = b >= inf;
    if inf1 && inf2 {
        return true;
    }
    if inf1 != inf2 {
        return false;
    }
    // assert( !inf1 ); assert( !inf2 );  // NDEBUG

    let minf1 = -a >= inf;
    let minf2 = -b >= inf;
    if minf1 && minf2 {
        return true;
    }
    if minf1 != minf2 {
        return false;
    }
    // assert( !minf1 ); assert( !minf2 );  // NDEBUG

    // EPSEQ(x,y,eps) = REALABS((x)-(y)) <= (eps)
    (a - b).abs() <= eps
}

/// SCIP 11.0.0 (32dd425): symmetry.c::SCIPsymLE (lines 1563-1597)
pub(crate) fn sym_le(a: f64, b: f64, eps: f64, inf: f64) -> bool {
    let inf1 = a >= inf;
    let inf2 = b >= inf;
    if inf1 && inf2 {
        return true;
    }
    if !inf1 && inf2 {
        return true;
    }
    if inf1 && !inf2 {
        return false;
    }
    // assert( !inf1 ); assert( !inf2 );  // NDEBUG

    let minf1 = -a >= inf;
    let minf2 = -b >= inf;
    if minf1 && minf2 {
        return true;
    }
    if !minf1 && minf2 {
        return false;
    }
    if minf1 && !minf2 {
        return true;
    }
    // assert( !minf1 ); assert( !minf2 );  // NDEBUG

    // EPSLE(x,y,eps) = (x)-(y) <= (eps)
    (a - b) <= eps
}

/// SCIP 11.0.0 (32dd425): symmetry.c::SCIPsymGE (lines 1601-1635)
pub(crate) fn sym_ge(a: f64, b: f64, eps: f64, inf: f64) -> bool {
    let inf1 = a >= inf;
    let inf2 = b >= inf;
    if inf1 && inf2 {
        return true;
    }
    if !inf1 && inf2 {
        return false;
    }
    if inf1 && !inf2 {
        return true;
    }
    // assert( !inf1 ); assert( !inf2 );  // NDEBUG

    let minf1 = -a >= inf;
    let minf2 = -b >= inf;
    if minf1 && minf2 {
        return true;
    }
    if !minf1 && minf2 {
        return true;
    }
    if minf1 && !minf2 {
        return false;
    }
    // assert( !minf1 ); assert( !minf2 );  // NDEBUG

    // EPSGE(x,y,eps) = (x)-(y) >= -(eps)
    (a - b) >= -eps
}

/// SCIP 11.0.0 (32dd425): symmetry.c::SCIPsymLT (lines 1639-1673)
pub(crate) fn sym_lt(a: f64, b: f64, eps: f64, inf: f64) -> bool {
    let inf1 = a >= inf;
    let inf2 = b >= inf;
    if inf1 && inf2 {
        return false;
    }
    if !inf1 && inf2 {
        return true;
    }
    if inf1 && !inf2 {
        return false;
    }
    // assert( !inf1 ); assert( !inf2 );  // NDEBUG

    let minf1 = -a >= inf;
    let minf2 = -b >= inf;
    if minf1 && minf2 {
        return false;
    }
    if !minf1 && minf2 {
        return false;
    }
    if minf1 && !minf2 {
        return true;
    }
    // assert( !minf1 ); assert( !minf2 );  // NDEBUG

    // EPSLT(x,y,eps) = (x)-(y) < -(eps)
    (a - b) < -eps
}

/// SCIP 11.0.0 (32dd425): symmetry.c::SCIPsymGT (lines 1677-1711)
pub(crate) fn sym_gt(a: f64, b: f64, eps: f64, inf: f64) -> bool {
    let inf1 = a >= inf;
    let inf2 = b >= inf;
    if inf1 && inf2 {
        return false;
    }
    if !inf1 && inf2 {
        return false;
    }
    if inf1 && !inf2 {
        return true;
    }
    // assert( !inf1 ); assert( !inf2 );  // NDEBUG

    let minf1 = -a >= inf;
    let minf2 = -b >= inf;
    if minf1 && minf2 {
        return false;
    }
    if !minf1 && minf2 {
        return true;
    }
    if minf1 && !minf2 {
        return false;
    }
    // assert( !minf1 ); assert( !minf2 );  // NDEBUG

    // EPSGT(x,y,eps) = (x)-(y) > (eps)
    (a - b) > eps
}

// ---------------------------------------------------------------------------
// bound-tightening wrappers
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): scip_var.c::SCIPtightenVarLb (lines 6447-6544)
///
/// Tighten a lower bound on `col` in `localdom`. Returns `(infeasible, tightened)`.
///
/// Same order as SCIP: round integer columns (`SCIPvarAdjustLb`: `ceil(val - feastol)`),
/// ignore a bound of +infinity, report infeasible if the new bound exceeds `ub + feastol`,
/// clamp to `ub`, skip a bound that does not tighten by more than `feastol`, and otherwise
/// call `change_bound_2` with an unspecified reason.
pub(crate) fn s3_tighten_lb(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    col: HighsInt,
    val: f64,
    is_int: bool,
    feastol: f64,
) -> (bool, bool) {
    let mut newval = if is_int { (val - feastol).ceil() } else { val };

    // ignore tightenings of lower bounds to +infinity during solving process
    if newval >= kHighsInf {
        return (false, false);
    }

    let col_lower = localdom.col_lower_[col as usize];
    let col_upper = localdom.col_upper_[col as usize];

    // SCIPsetIsFeasGT(newbound, ub)
    if newval > col_upper + feastol {
        return (true, false);
    }
    newval = cpp_min(newval, col_upper);

    // !SCIPsetIsLbBetter(newbound, lb, ub), D2: tighten by more than feastol
    if newval <= col_lower + feastol {
        return (false, false);
    }

    localdom.change_bound_2(
        ms,
        HighsBoundType::kLower,
        col,
        newval,
        HighsDomain_Reason::unspecified(),
    );

    (localdom.infeasible(), true)
}

/// SCIP 11.0.0 (32dd425): scip_var.c::SCIPtightenVarUb (lines 6697-6794)
///
/// Tighten an upper bound on `col` in `localdom`. Returns `(infeasible, tightened)`.
///
/// Mirror of `s3_tighten_lb`: round integer columns (`floor(val + feastol)`), ignore a
/// bound of -infinity, report infeasible if the new bound is below `lb - feastol`, clamp
/// to `lb`, skip a bound that does not tighten by more than `feastol`, and otherwise call
/// `change_bound_2` with an unspecified reason.
pub(crate) fn s3_tighten_ub(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    col: HighsInt,
    val: f64,
    is_int: bool,
    feastol: f64,
) -> (bool, bool) {
    let mut newval = if is_int { (val + feastol).floor() } else { val };

    // ignore tightenings of upper bounds to -infinity during solving process
    if -newval >= kHighsInf {
        return (false, false);
    }

    let col_lower = localdom.col_lower_[col as usize];
    let col_upper = localdom.col_upper_[col as usize];

    // SCIPsetIsFeasLT(newbound, lb)
    if newval < col_lower - feastol {
        return (true, false);
    }
    newval = cpp_max(newval, col_lower);

    // !SCIPsetIsUbBetter(newbound, lb, ub), D2: tighten by more than feastol
    if newval >= col_upper - feastol {
        return (false, false);
    }

    localdom.change_bound_2(
        ms,
        HighsBoundType::kUpper,
        col,
        newval,
        HighsDomain_Reason::unspecified(),
    );

    (localdom.infeasible(), true)
}

// ---------------------------------------------------------------------------
// branching path
// ---------------------------------------------------------------------------

/// Per-orbitope information recorded at a branching node.
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c (nodeinfo entries, lines 940-960)
#[derive(Debug, Clone, Default)]
pub(crate) struct S3OrbiNodeInfo {
    pub(crate) orbitope: u32,
    pub(crate) rows_added: Vec<HighsInt>,
    pub(crate) colswaps: Vec<(HighsInt, HighsInt)>,
}

/// A node in the branching path, replacing SCIP's nodeinfo hash table and
/// shadow tree. Both children of a branching share the same `Arc<S3PathNode>`.
///
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::eventExecNodeBranched (lines 720-962)
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c (shadow tree, lines 1640-1850)
#[derive(Debug)]
pub(crate) struct S3PathNode {
    pub(crate) parent: Option<Arc<S3PathNode>>,
    pub(crate) branch_col: HighsInt,
    pub(crate) depth: HighsInt,
    pub(crate) orbi: Vec<S3OrbiNodeInfo>,
}

// ---------------------------------------------------------------------------
// S3OrbitopeData
// ---------------------------------------------------------------------------

/// Orbitopal symmetry handling data for a single orbitope.
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::OrbitopeData (lines 126-142)
///
/// Fields `nodeinfos`, `lastnodenumber` and `dbghash` are NOT ported (R5, D3).
#[derive(Debug, Clone, Default)]
pub(crate) struct S3OrbitopeData {
    /// HiGHS column indices, indexed [row][col].
    pub(crate) vars: Vec<Vec<HighsInt>>,
    pub(crate) nrows: HighsInt,
    pub(crate) ncols: HighsInt,
    pub(crate) nbranchrows: HighsInt,
    /// Map from variable index to row index in the orbitope matrix.
    pub(crate) rowindexmap: Vec<HighsInt>,
    /// Map from variable index to column index in the orbitope matrix.
    pub(crate) colindexmap: Vec<HighsInt>,
    pub(crate) columnordering: HighsInt,
    pub(crate) rowordering: HighsInt,
}

// ---------------------------------------------------------------------------
// S3LexData
// ---------------------------------------------------------------------------

/// Lexicographic reduction permutation data.
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::LexRedPermData (lines 89-100)
#[derive(Debug, Clone, Default)]
pub(crate) struct S3LexData {
    pub(crate) isdynamic: bool,
    pub(crate) vars: Vec<HighsInt>,
    pub(crate) perm: Vec<HighsInt>,
    pub(crate) invperm: Vec<HighsInt>,
    /// Lookup only, never iterated.
    pub(crate) varmap: std::collections::HashMap<HighsInt, HighsInt>,
    pub(crate) symvarmap: Vec<HighsInt>,
    pub(crate) nsymvars: HighsInt,
    pub(crate) hasdynamicperm: bool,
}

// ---------------------------------------------------------------------------
// S3SymData
// ---------------------------------------------------------------------------

/// Top-level container for all symmetry-handling data.
#[derive(Debug, Clone, Default)]
pub(crate) struct S3SymData {
    pub(crate) orbitopes: Vec<S3OrbitopeData>,
    pub(crate) static_orbitopes: Vec<usize>,
    pub(crate) lexdatas: Vec<S3LexData>,
}

// ---------------------------------------------------------------------------
// S3Stats
// ---------------------------------------------------------------------------

/// Per-propagator statistics counters.
#[derive(Debug, Clone, Default)]
pub(crate) struct S3Stats {
    pub(crate) orbitope_nred: u64,
    pub(crate) orbitope_ncutoff: u64,
    pub(crate) lexred_nred: u64,
    pub(crate) lexred_ncutoff: u64,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s3_data_sym_comparisons_with_infinities() {
        let eps = 1e-9;
        let inf = kHighsInf;

        // Two positive infinities are equal
        assert!(sym_eq(f64::INFINITY, f64::INFINITY, eps, inf));
        // Two negative infinities are equal
        assert!(sym_eq(f64::NEG_INFINITY, f64::NEG_INFINITY, eps, inf));
        // +inf != -inf
        assert!(!sym_eq(f64::INFINITY, f64::NEG_INFINITY, eps, inf));
        // inf != finite
        assert!(!sym_eq(f64::INFINITY, 1.0, eps, inf));
        assert!(!sym_eq(1.0, f64::INFINITY, eps, inf));

        // LE: inf <= inf is true
        assert!(sym_le(f64::INFINITY, f64::INFINITY, eps, inf));
        // LE: finite <= inf is true
        assert!(sym_le(1.0, f64::INFINITY, eps, inf));
        // LE: inf <= finite is false
        assert!(!sym_le(f64::INFINITY, 1.0, eps, inf));
        // LE: -inf <= -inf is true
        assert!(sym_le(f64::NEG_INFINITY, f64::NEG_INFINITY, eps, inf));
        // LE: -inf <= finite is true
        assert!(sym_le(f64::NEG_INFINITY, 1.0, eps, inf));
        // LE: finite <= -inf is false
        assert!(!sym_le(1.0, f64::NEG_INFINITY, eps, inf));

        // GE: inf >= inf is true
        assert!(sym_ge(f64::INFINITY, f64::INFINITY, eps, inf));
        // GE: inf >= finite is true
        assert!(sym_ge(f64::INFINITY, 1.0, eps, inf));
        // GE: finite >= inf is false
        assert!(!sym_ge(1.0, f64::INFINITY, eps, inf));
        // GE: -inf >= -inf is true
        assert!(sym_ge(f64::NEG_INFINITY, f64::NEG_INFINITY, eps, inf));
        // GE: -inf >= finite is false
        assert!(!sym_ge(f64::NEG_INFINITY, 1.0, eps, inf));
        // GE: finite >= -inf is true
        assert!(sym_ge(1.0, f64::NEG_INFINITY, eps, inf));

        // LT: inf < inf is false
        assert!(!sym_lt(f64::INFINITY, f64::INFINITY, eps, inf));
        // LT: finite < inf is true
        assert!(sym_lt(1.0, f64::INFINITY, eps, inf));
        // LT: inf < finite is false
        assert!(!sym_lt(f64::INFINITY, 1.0, eps, inf));
        // LT: -inf < -inf is false
        assert!(!sym_lt(f64::NEG_INFINITY, f64::NEG_INFINITY, eps, inf));
        // LT: -inf < finite is false (minf1 && !minf2 -> true? Let's check C code)
        // In SCIPsymLT: minf1 && !minf2 -> return TRUE. So -inf < finite is TRUE.
        assert!(sym_lt(f64::NEG_INFINITY, 1.0, eps, inf));
        // LT: finite < -inf is false
        assert!(!sym_lt(1.0, f64::NEG_INFINITY, eps, inf));

        // GT: inf > inf is false
        assert!(!sym_gt(f64::INFINITY, f64::INFINITY, eps, inf));
        // GT: inf > finite is true
        assert!(sym_gt(f64::INFINITY, 1.0, eps, inf));
        // GT: finite > inf is false
        assert!(!sym_gt(1.0, f64::INFINITY, eps, inf));
        // GT: -inf > -inf is false
        assert!(!sym_gt(f64::NEG_INFINITY, f64::NEG_INFINITY, eps, inf));
        // GT: -inf > finite is false
        assert!(!sym_gt(f64::NEG_INFINITY, 1.0, eps, inf));
        // GT: finite > -inf is true
        assert!(sym_gt(1.0, f64::NEG_INFINITY, eps, inf));
    }

    #[test]
    fn s3_data_sym_comparisons_with_epsilon() {
        let eps = 1e-9;
        let inf = kHighsInf;

        // EQ: values within eps are equal
        assert!(sym_eq(1.0, 1.0 + eps / 2.0, eps, inf));
        assert!(sym_eq(1.0, 1.0 - eps / 2.0, eps, inf));
        // EQ: values outside eps are not equal
        assert!(!sym_eq(1.0, 1.0 + eps * 2.0, eps, inf));

        // LE: a <= b + eps
        assert!(sym_le(1.0, 1.0, eps, inf));
        assert!(sym_le(1.0, 1.0 + eps / 2.0, eps, inf));
        assert!(!sym_le(1.0 + eps * 2.0, 1.0, eps, inf));

        // GE: a >= b - eps
        assert!(sym_ge(1.0, 1.0, eps, inf));
        assert!(sym_ge(1.0, 1.0 - eps / 2.0, eps, inf));
        assert!(!sym_ge(1.0 - eps * 2.0, 1.0, eps, inf));

        // LT: a < b - eps
        assert!(!sym_lt(1.0, 1.0, eps, inf));
        assert!(sym_lt(1.0 - eps * 2.0, 1.0, eps, inf));
        assert!(!sym_lt(1.0, 1.0 + eps / 2.0, eps, inf));

        // GT: a > b + eps
        assert!(!sym_gt(1.0, 1.0, eps, inf));
        assert!(sym_gt(1.0 + eps * 2.0, 1.0, eps, inf));
        assert!(!sym_gt(1.0, 1.0 - eps / 2.0, eps, inf));
    }

    /// A one-column integer MIP (no rows), so that `change_bound` finds a model behind `ms`.
    fn one_col_solver() -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 1;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![10.0];
        lp.integrality_ = vec![HighsVarType::kInteger];
        lp.a_matrix_.num_col_ = 1;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0, 0];
        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0)
    }

    #[test]
    fn s3_data_tighten_lb_rounds_integer_bounds() {
        let feastol = 1e-6;
        let mut ms = one_col_solver();
        let mut localdom = HighsDomain::default();
        localdom.col_lower_.resize(1, 0.0);
        localdom.col_upper_.resize(1, 10.0);
        localdom.col_lower_pos_ = vec![-1; 1];
        localdom.col_upper_pos_ = vec![-1; 1];
        localdom.changedcolsflags_ = vec![0; 1];

        // Integer column: val=3.7, feastol=1e-6
        // newval = ceil(3.7 - 1e-6) = ceil(3.699999) = 4
        // 4 > 0 + 1e-6, so it tightens
        let (infeasible, tightened) = s3_tighten_lb(&mut ms, &mut localdom, 0, 3.7, true, feastol);
        assert!(!infeasible);
        assert!(tightened);
        assert_eq!(localdom.col_lower_[0], 4.0);
    }

    #[test]
    fn s3_data_tighten_ub_skips_non_tightening() {
        let feastol = 1e-6;
        let mut ms = HighsMipSolver::default();
        let mut localdom = HighsDomain::default();
        localdom.col_lower_.resize(1, 0.0);
        localdom.col_upper_.resize(1, 10.0);
        localdom.col_lower_pos_ = vec![-1; 1];
        localdom.col_upper_pos_ = vec![-1; 1];
        localdom.changedcolsflags_ = vec![0; 1];

        // Upper bound: val=10.0, feastol=1e-6
        // newval = 10.0 (continuous)
        // 10.0 >= 10.0 - 1e-6, so it does NOT tighten
        let (infeasible, tightened) =
            s3_tighten_ub(&mut ms, &mut localdom, 0, 10.0, false, feastol);
        assert!(!infeasible);
        assert!(!tightened);
        assert_eq!(localdom.col_upper_[0], 10.0);

        // SCIP ignores tightenings of an upper bound to -infinity and of a lower bound
        // to +infinity.
        let (infeasible, tightened) =
            s3_tighten_ub(&mut ms, &mut localdom, 0, -kHighsInf, true, feastol);
        assert!(!infeasible);
        assert!(!tightened);
        let (infeasible, tightened) =
            s3_tighten_lb(&mut ms, &mut localdom, 0, kHighsInf, true, feastol);
        assert!(!infeasible);
        assert!(!tightened);
    }

    #[test]
    fn s3_data_tighten_detects_infeasible() {
        let feastol = 1e-6;
        let mut ms = HighsMipSolver::default();
        let mut localdom = HighsDomain::default();
        localdom.col_lower_.resize(1, 5.0);
        localdom.col_upper_.resize(1, 3.0);
        localdom.col_lower_pos_ = vec![-1; 1];
        localdom.col_upper_pos_ = vec![-1; 1];
        localdom.changedcolsflags_ = vec![0; 1];

        // Lower bound: val=6.0, feastol=1e-6
        // newval = 6.0 (continuous)
        // 6.0 > 3.0 + 1e-6, so infeasible
        let (infeasible, tightened) = s3_tighten_lb(&mut ms, &mut localdom, 0, 6.0, false, feastol);
        assert!(infeasible);
        assert!(!tightened);
    }

    #[test]
    fn s3_data_s3_path_node_shares_parent() {
        let root = Arc::new(S3PathNode {
            parent: None,
            branch_col: -1,
            depth: 0,
            orbi: vec![],
        });

        let child = Arc::new(S3PathNode {
            parent: Some(root.clone()),
            branch_col: 5,
            depth: 1,
            orbi: vec![],
        });

        // Both children share the same parent Arc
        let child2 = Arc::new(S3PathNode {
            parent: Some(root.clone()),
            branch_col: 7,
            depth: 1,
            orbi: vec![],
        });

        assert!(Arc::ptr_eq(
            child.parent.as_ref().unwrap(),
            child2.parent.as_ref().unwrap()
        ));
        assert_eq!(child.branch_col, 5);
        assert_eq!(child2.branch_col, 7);
        assert_eq!(child.depth, 1);
    }
}
