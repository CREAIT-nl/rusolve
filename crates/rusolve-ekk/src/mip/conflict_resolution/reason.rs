//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c (isResolvableBdchg 1012-1037, reasonRowFromLpRow 2306-2396,
//!   conflictRowFromLpRow 2398, getReasonRow 2481-2568); HiGHS HighsDomain.cpp (explainInfeasibility
//!   3205-3302, explainInfeasibilityLeq 3423-3467, conflicting bounds 2040-2070, getCol*Pos 2563-2583).

use super::*;
use crate::prelude::*;

/// S2Context implementation over a HiGHS MIP solver and its domains.
pub struct S2DomainCtx<'a> {
    pub ms: &'a HighsMipSolver,
    pub localdom: &'a mut HighsDomain,
    pub globaldom: &'a HighsDomain,
}

impl<'a> S2Context for S2DomainCtx<'a> {
    fn num_col(&self) -> usize {
        self.localdom.col_lower_.len()
    }

    fn num_integral(&self) -> usize {
        self.ms.mipdata_.integral_cols.len()
    }

    fn is_integral(&self, col: HighsInt) -> bool {
        self.ms.is_col_integral(col)
    }

    fn global_lb(&self, col: HighsInt) -> f64 {
        self.globaldom.col_lower_[col as usize]
    }

    fn global_ub(&self, col: HighsInt) -> f64 {
        self.globaldom.col_upper_[col as usize]
    }

    fn lb_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt) {
        let mut pos = -1;
        let bound = self.localdom.get_col_lower_pos(col, upto, &mut pos);
        (bound, pos)
    }

    fn ub_upto(&self, col: HighsInt, upto: HighsInt) -> (f64, HighsInt) {
        let mut pos = -1;
        let bound = self.localdom.get_col_upper_pos(col, upto, &mut pos);
        (bound, pos)
    }

    fn change(&self, pos: HighsInt) -> S2Change {
        let chg = &self.localdom.domchgstack_[pos as usize];
        let reason = self.localdom.domchgreason_[pos as usize];
        let col = chg.column;
        let is_upper = chg.boundtype == HighsBoundType::kUpper;
        let newbound = chg.boundval;
        let oldbound = self.localdom.prevboundval_[pos as usize].0;
        let depth = self.localdom.branch_pos_.partition_point(|&b| b <= pos) as HighsInt;
        let branching = reason.r#type == HighsDomain_Reason::kBranching;
        let resolvable = s2_reason_is_linear(self, pos);
        S2Change {
            col,
            is_upper,
            newbound,
            oldbound,
            depth,
            branching,
            resolvable,
        }
    }

    fn reason_row(&mut self, pos: HighsInt, out: &mut S2ConflictRow) -> bool {
        let chg = &self.localdom.domchgstack_[pos as usize];
        let col = chg.column;
        let is_upper = chg.boundtype == HighsBoundType::kUpper;
        let reason = self.localdom.domchgreason_[pos as usize];
        self.materialise(reason, pos, col, is_upper, out)
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::isResolvableBdchg (lines 1012-1037)
/// HiGHS: HighsDomain.cpp::ConflictSet::explainInfeasibility (lines 3205-3302)
pub fn s2_reason_is_linear(ctx: &S2DomainCtx<'_>, pos: HighsInt) -> bool {
    let reason = ctx.localdom.domchgreason_[pos as usize];
    match reason.r#type {
        HighsDomain_Reason::kBranching
        | HighsDomain_Reason::kUnknown
        | HighsDomain_Reason::kConflictingBounds => false,
        HighsDomain_Reason::kModelRowUpper
        | HighsDomain_Reason::kModelRowLower
        | HighsDomain_Reason::kObjective => true,
        HighsDomain_Reason::kCliqueTable => {
            let chg = &ctx.localdom.domchgstack_[pos as usize];
            s2_clique_change_is_fixing(
                ctx,
                chg.column,
                chg.boundtype == HighsBoundType::kUpper,
                chg.boundval,
            )
        }
        t if t >= 0 => {
            let ncutpools = ctx.localdom.cutpoolpropagation.len() as HighsInt;
            if t < ncutpools {
                true
            } else {
                let cpidx = (t - ncutpools) as usize;
                debug_assert!(cpidx < ctx.localdom.conflict_pool_propagation.len());
                let cp = &ctx.localdom.conflict_pool_propagation[cpidx];
                let pool = &ctx.ms.mipdata_.conflictpools[cp.conflictpool];
                if cp.conflict_flag_[reason.index as usize] & 8 != 0 {
                    return false;
                }
                let ranges = pool.get_conflict_ranges();
                let entries = pool.get_conflict_entry_vector();
                let (start, end) = ranges[reason.index as usize];
                for e in &entries[start as usize..end as usize] {
                    if !s2_is_binary(ctx, e.column) {
                        return false;
                    }
                    if e.boundtype == HighsBoundType::kLower && e.boundval != 1.0 {
                        return false;
                    }
                    if e.boundtype == HighsBoundType::kUpper && e.boundval != 0.0 {
                        return false;
                    }
                }
                true
            }
        }
        _ => false,
    }
}

impl<'a> S2DomainCtx<'a> {
    /// HiGHS: HighsDomain.cpp::ConflictSet::explainInfeasibility (lines 3205-3302)
    pub fn materialise(
        &mut self,
        reason: HighsDomain_Reason,
        stackpos: HighsInt,
        col: HighsInt,
        is_upper: bool,
        out: &mut S2ConflictRow,
    ) -> bool {
        let ncol = self.num_col();
        out.clear(ncol);

        let ok = match reason.r#type {
            HighsDomain_Reason::kModelRowUpper => {
                let r = reason.index;
                let mut len = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                self.ms.get_row(r, &mut len, &mut inds, &mut vals);
                for k in 0..len as usize {
                    out.vals[inds[k] as usize] = -vals[k];
                    out.inds.push(inds[k]);
                }
                out.lhs = -self.ms.row_upper_2(r);
                true
            }
            HighsDomain_Reason::kModelRowLower => {
                let r = reason.index;
                let mut len = 0;
                let mut inds: &[HighsInt] = &[];
                let mut vals: &[f64] = &[];
                self.ms.get_row(r, &mut len, &mut inds, &mut vals);
                for k in 0..len as usize {
                    out.vals[inds[k] as usize] = vals[k];
                    out.inds.push(inds[k]);
                }
                out.lhs = self.ms.row_lower_2(r);
                true
            }
            HighsDomain_Reason::kCliqueTable => {
                // the D3 clique rows hold only for a binary fixing of y (see
                // s2_clique_change_is_fixing); implication reasons on other bounds are not linear
                let boundval = self.localdom.domchgstack_[stackpos as usize].boundval;
                if !s2_clique_change_is_fixing(&*self, col, is_upper, boundval) {
                    return false;
                }
                let xcol = reason.index >> 1;
                let val = reason.index & 1;
                clique_reason_row(xcol, val, col, is_upper, out);
                true
            }
            HighsDomain_Reason::kObjective => {
                let mut vals: &[f64] = &[];
                let mut inds: &[HighsInt] = &[];
                let mut len = 0;
                let mut rhs = 0.0;
                self.localdom.obj_prop_get_propagation_constraint(
                    self.ms, stackpos, &mut vals, &mut inds, &mut len, &mut rhs, col,
                );
                let vals_vec = vals.to_vec();
                let inds_vec = inds.to_vec();
                // the row is built directly; rows with infinite bounds are rejected by
                // the slack tests (getReasonRow / getConflictRow, ), not here.
                for k in 0..len as usize {
                    out.vals[inds_vec[k] as usize] = -vals_vec[k];
                    out.inds.push(inds_vec[k]);
                }
                out.lhs = -rhs;
                true
            }
            t if t >= 0 => {
                let ncutpools = self.localdom.cutpoolpropagation.len() as HighsInt;
                if t < ncutpools {
                    let pool_idx = self.localdom.cutpoolpropagation[t as usize].cutpool;
                    let pool = &self.ms.mipdata_.cutpools[pool_idx];
                    let mut len = 0;
                    let mut inds: &[HighsInt] = &[];
                    let mut vals: &[f64] = &[];
                    pool.get_cut(reason.index, &mut len, &mut inds, &mut vals);
                    for k in 0..len as usize {
                        out.vals[inds[k] as usize] = -vals[k];
                        out.inds.push(inds[k]);
                    }
                    out.lhs = -pool.get_rhs()[reason.index as usize];
                    true
                } else {
                    let cpidx = (t - ncutpools) as usize;
                    debug_assert!(cpidx < self.localdom.conflict_pool_propagation.len());
                    let cp = &self.localdom.conflict_pool_propagation[cpidx];
                    let pool = &self.ms.mipdata_.conflictpools[cp.conflictpool];
                    if cp.conflict_flag_[reason.index as usize] & 8 != 0 {
                        return false;
                    }
                    let ranges = pool.get_conflict_ranges();
                    let entries = pool.get_conflict_entry_vector();
                    let (start, end) = ranges[reason.index as usize];
                    let mut nogood_entries = Vec::new();
                    for e in &entries[start as usize..end as usize] {
                        let is_up = e.boundtype == HighsBoundType::kUpper;
                        nogood_entries.push((e.column, is_up, e.boundval));
                    }
                    let is_bin = |col: HighsInt| s2_is_binary(&*self, col);
                    nogood_reason_row(&nogood_entries, is_bin, out)
                }
            }
            _ => false,
        };
        // Rows have unique indices (manifest the port).
        debug_assert!(!ok || s2_inds_unique(&out.inds));
        ok && s2_drop_zero_coefs(&*self, out)
    }
}

/// SCIP conflict rows never hold a zero coefficient. SCIP builds them from LP rows and
/// linear constraints, which store none (reasonRowFromLpRow 2306-2396), and linearCombRows
/// (975-1008) relies on that: it appends `idx` whenever `vals[idx]` is zero, so a column kept
/// in `inds` with a zero coefficient gets appended a second time. HiGHS rows can hold zeros:
/// `ObjectivePropagation::getPropagationConstraint` (HighsDomain.cpp:780-786) sets a clique
/// partition's coefficients to `packedObjVals[j] - copysign(largest, packedObjVals[j])`, which
/// is exactly 0 for every column whose cost equals the partition's largest cost.
///
/// Exact zeros are dropped, which leaves the row unchanged. A coefficient within S2_EPSILON of
/// zero is weakened out with its global bound (`a > 0`: lhs -= a*ub, `a < 0`: lhs -= a*lb). This
/// keeps the `>=` row valid. If that bound is infinite the row cannot be used.
pub fn s2_drop_zero_coefs<C: S2Context>(ctx: &C, row: &mut S2ConflictRow) -> bool {
    for i in (0..row.nnz()).rev() {
        let idx = row.inds[i];
        let v = row.vals[idx as usize];
        if v != 0.0 && s2_is_zero(v) {
            let bnd = if v > 0.0 {
                ctx.global_ub(idx)
            } else {
                ctx.global_lb(idx)
            };
            if s2_is_infinity(bnd.abs()) {
                return false;
            }
            row.lhs -= v * bnd;
            remove_zero_var(row, i);
        } else if v == 0.0 {
            remove_zero_var(row, i);
        }
    }
    true
}

/// Debug helper: every column appears at most once in `inds`.
fn s2_inds_unique(inds: &[HighsInt]) -> bool {
    let mut sorted = inds.to_vec();
    sorted.sort_unstable();
    sorted.windows(2).all(|w| w[0] != w[1])
}

/// a `kCliqueTable` reason gives the clique row `lit(x,val) + lit(y) <= 1` only when the
/// change fixes a binary y (ub to 0 or lb to 1). HiGHS also uses `kCliqueTable` for probing
/// implications (`$HIGHS`/mip/HighsImplications.cpp:877-905, `applyImplications`), which can
/// tighten a general-integer or continuous y to any value; the clique row is invalid for those,
/// so such changes are not resolvable (SCIP: implication changes are not CONSINFER,
/// conflict_resolution.c:1012-1037).
pub fn s2_clique_change_is_fixing<C: S2Context>(
    ctx: &C,
    ycol: HighsInt,
    is_upper: bool,
    boundval: f64,
) -> bool {
    s2_is_binary(ctx, ycol)
        && if is_upper {
            boundval == 0.0
        } else {
            boundval == 1.0
        }
}

/// clique reason row.
pub fn clique_reason_row(
    xcol: HighsInt,
    val: HighsInt,
    ycol: HighsInt,
    is_upper: bool,
    out: &mut S2ConflictRow,
) {
    // Clear `out`, keeping it dense over at least both columns.
    let ncol = cpp_max(out.vals.len(), (cpp_max(xcol, ycol) + 1) as usize);
    out.clear(ncol);

    if is_upper {
        if val != 0 {
            out.vals[xcol as usize] = -1.0;
            out.vals[ycol as usize] = -1.0;
            out.lhs = -1.0;
        } else {
            out.vals[xcol as usize] = 1.0;
            out.vals[ycol as usize] = -1.0;
            out.lhs = 0.0;
        }
    } else {
        if val != 0 {
            out.vals[xcol as usize] = -1.0;
            out.vals[ycol as usize] = 1.0;
            out.lhs = 0.0;
        } else {
            out.vals[xcol as usize] = 1.0;
            out.vals[ycol as usize] = 1.0;
            out.lhs = 1.0;
        }
    }
    out.inds.push(xcol);
    out.inds.push(ycol);
}

/// nogood reason row from conflict pool entries.
pub fn nogood_reason_row(
    entries: &[(HighsInt, bool, f64)],
    is_binary: impl Fn(HighsInt) -> bool,
    out: &mut S2ConflictRow,
) -> bool {
    for &(col, is_upper, boundval) in entries {
        if !is_binary(col) {
            return false;
        }
        if !is_upper && boundval != 1.0 {
            return false;
        }
        if is_upper && boundval != 0.0 {
            return false;
        }
    }

    let mut max_col: HighsInt = -1;
    for &(col, _, _) in entries {
        if col > max_col {
            max_col = col;
        }
    }
    // Clear `out`, keeping it dense over at least every entry column.
    let ncol = cpp_max(out.vals.len(), (max_col + 1) as usize);
    out.clear(ncol);

    let mut num_lower = 0;
    for &(col, is_upper, _) in entries {
        if is_upper {
            out.vals[col as usize] = 1.0;
        } else {
            out.vals[col as usize] = -1.0;
            num_lower += 1;
        }
        out.inds.push(col);
    }
    out.lhs = 1.0 - num_lower as f64;
    true
}

/// initial conflict row from proof or infeasibility reason.
pub fn s2_initial_conflict_row(
    ctx: &mut S2DomainCtx<'_>,
    proof: Option<(&[HighsInt], &[f64], HighsInt, f64)>,
    out: &mut S2ConflictRow,
) -> Option<HighsInt> {
    if let Some((inds, vals, len, rhs)) = proof {
        let ncol = ctx.num_col();
        out.clear(ncol);
        for k in 0..len as usize {
            out.vals[inds[k] as usize] = -vals[k];
            out.inds.push(inds[k]);
        }
        out.lhs = -rhs;
        if !s2_drop_zero_coefs(&*ctx, out) {
            return None;
        }
        return Some(kHighsIInf);
    }

    // infeasible_reason is only defined when the domain is infeasible; HighsDomain::conflictAnalysis
    // returns on `!infeasible_` before analysis (stack.rs conflict_analysis).
    if !ctx.localdom.infeasible_ {
        return None;
    }

    let reason = ctx.localdom.infeasible_reason;
    match reason.r#type {
        HighsDomain_Reason::kUnknown
        | HighsDomain_Reason::kBranching
        | HighsDomain_Reason::kCliqueTable => None,
        HighsDomain_Reason::kConflictingBounds => {
            let p1 = reason.index;
            if ctx.reason_row(p1, out) {
                return Some(p1 - 1);
            }
            // HighsDomain.cpp:3214-3226: p2 is the latest change of the opposite bound.
            let chg1 = &ctx.localdom.domchgstack_[p1 as usize];
            let col = chg1.column;
            let mut p2: HighsInt = -1;
            if chg1.boundtype == HighsBoundType::kLower {
                ctx.localdom.get_col_upper_pos(col, p1, &mut p2);
            } else {
                ctx.localdom.get_col_lower_pos(col, p1, &mut p2);
            }
            if p2 != -1 && ctx.reason_row(p2, out) {
                Some(p1)
            } else {
                None
            }
        }
        _ => {
            let pos = ctx.localdom.infeasible_pos;
            if ctx.materialise(reason, pos, -1, false, out) {
                Some(pos)
            } else {
                None
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::getReasonRow (lines 2481-2568)
pub fn get_reason_row(
    ctx: &mut impl S2Context,
    pos: HighsInt,
    reasonrow: &mut S2ConflictRow,
    residx: HighsInt,
    fixbounds: &[f64],
    fixsides: &[i8],
) -> bool {
    let change = ctx.change(pos);
    if !change.resolvable {
        return false;
    }
    if !ctx.reason_row(pos, reasonrow) {
        return false;
    }
    // SCIP knapsack negated-clique test (conflict_resolution.c 2527-2532) is dropped:
    // HiGHS reason types already name the side.
    compute_slack(&*ctx, reasonrow, pos, Some((fixbounds, fixsides)));
    if s2_is_gt(reasonrow.slack, 0.0) {
        let coef = reasonrow.vals[residx as usize];
        let boundusedinslack = if coef > 0.0 {
            s2_ub_at(&*ctx, change.col, pos, true)
        } else {
            s2_lb_at(&*ctx, change.col, pos, true)
        };
        if !s2_is_lt(
            reasonrow.slack - coef * (boundusedinslack - change.oldbound),
            0.0,
        ) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_mip() -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 2;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![0.0, 0.0];
        lp.col_lower_ = vec![0.0, 0.0];
        lp.col_upper_ = vec![10.0, 10.0];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![4.0];
        lp.integrality_ = vec![HighsVarType::kInteger; 2];
        lp.a_matrix_.num_col_ = 2;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2];
        lp.a_matrix_.index_ = vec![0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0];
        ms.mipdata_.a_rstart_ = vec![0, 2];
        ms.mipdata_.a_rindex_ = vec![0, 1];
        ms.mipdata_.a_rvalue_ = vec![1.0, 1.0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms
    }

    fn lower(col: HighsInt, val: f64) -> HighsDomainChange {
        HighsDomainChange {
            boundval: val,
            column: col,
            boundtype: HighsBoundType::kLower,
        }
    }

    #[test]
    fn s2_reason_model_row() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);
        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        assert!(dom.propagate(&mut ms));
        assert!(!dom.infeasible());
        assert_eq!(dom.col_upper_[1], 1.0);
        // Propagation records x0 <= 4 at pos 1 and x1 <= 1 at pos 2 (both kModelRowUpper(0));
        // the manifest's "pos 1" counted only the x1 change.
        assert_eq!(dom.domchgstack_[2].column, 1);
        assert_eq!(dom.domchgstack_[2].boundtype, HighsBoundType::kUpper);
        assert_eq!(
            dom.domchgreason_[2].r#type,
            HighsDomain_Reason::kModelRowUpper
        );

        let globaldom = HighsDomain::new(&ms);
        let mut ctx = S2DomainCtx {
            ms: &ms,
            localdom: &mut dom,
            globaldom: &globaldom,
        };

        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        assert!(ctx.reason_row(2, &mut row));
        assert_eq!(row.vals, vec![-1.0, -1.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, -4.0);

        let fixbounds = vec![0.0; 2];
        let fixsides = vec![0i8; 2];
        assert!(get_reason_row(
            &mut ctx, 2, &mut row, 1, &fixbounds, &fixsides
        ));
        assert_eq!(row.slack, 1.0);
    }

    #[test]
    fn s2_reason_conflicting_bounds() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);
        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        assert!(dom.propagate(&mut ms));
        dom.change_bound(&mut ms, lower(1, 2.0), HighsDomain_Reason::branching());
        assert!(dom.infeasible_);
        assert_eq!(
            dom.infeasible_reason.r#type,
            HighsDomain_Reason::kConflictingBounds
        );
        // Propagating x0 >= 3 also tightens x0 <= 4 besides x1 <= 1, so the stack is
        // [x0>=3 (branch), two kModelRowUpper(0) changes, x1>=2 (branch)] and p1 = 3
        // (the manifest's p1 = 2 counted only the x1 change).
        let p1 = dom.infeasible_reason.index;
        assert_eq!(p1, 3);
        let mut p2: HighsInt = -1;
        dom.get_col_upper_pos(1, p1, &mut p2);
        assert_eq!(p2, 2);
        assert_eq!(
            dom.domchgreason_[p2 as usize].r#type,
            HighsDomain_Reason::kModelRowUpper
        );

        let globaldom = HighsDomain::new(&ms);
        let mut ctx = S2DomainCtx {
            ms: &ms,
            localdom: &mut dom,
            globaldom: &globaldom,
        };

        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        // reason(p1) is a branching, so case (b): R2 = row of reason(p2), upto = p1.
        let upto = s2_initial_conflict_row(&mut ctx, None, &mut row);
        assert_eq!(upto, Some(p1));
        assert_eq!(row.vals, vec![-1.0, -1.0]);
        assert_eq!(row.lhs, -4.0);
        compute_slack(&ctx, &mut row, p1, None);
        assert_eq!(row.slack, -1.0);
    }

    #[test]
    fn s2_reason_row_infeasible() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        dom.compute_row_activities(&ms);
        dom.change_bound(&mut ms, lower(0, 3.0), HighsDomain_Reason::branching());
        dom.change_bound(&mut ms, lower(1, 2.0), HighsDomain_Reason::branching());
        assert!(dom.propagate(&mut ms));
        assert!(dom.infeasible_);
        assert_eq!(
            dom.infeasible_reason.r#type,
            HighsDomain_Reason::kModelRowUpper
        );
        assert_eq!(dom.infeasible_reason.index, 0);
        // HighsDomain.cpp:1576-1585 sets infeasible_pos = domchgstack_.size() inside
        // doChangeBound, before the push: the activity check fires while x1 >= 2 is being
        // applied, so infeasible_pos is that change's position 1 (the manifest's "2"
        // assumed a later detection). getColLowerPos(col, 1) still counts both branchings.
        assert_eq!(dom.infeasible_pos, 1);

        let globaldom = HighsDomain::new(&ms);
        let mut ctx = S2DomainCtx {
            ms: &ms,
            localdom: &mut dom,
            globaldom: &globaldom,
        };

        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        let upto = s2_initial_conflict_row(&mut ctx, None, &mut row);
        assert_eq!(upto, Some(1));
        assert_eq!(row.vals, vec![-1.0, -1.0]);
        assert_eq!(row.lhs, -4.0);
        compute_slack(&ctx, &mut row, 1, None);
        assert_eq!(row.slack, -1.0);
    }

    #[test]
    fn s2_reason_clique_implication_not_fixing() {
        // col 0 binary, col 1 general integer [0, 5], col 2 continuous [0, 1]
        let ctx = crate::mip::conflict_resolution::S2TestCtx {
            glb: vec![0.0, 0.0, 0.0],
            gub: vec![1.0, 5.0, 1.0],
            integral: vec![true, true, false],
            num_integral: 2,
            stack: vec![],
        };
        // binary fixings: the clique row holds
        assert!(s2_clique_change_is_fixing(&ctx, 0, true, 0.0));
        assert!(s2_clique_change_is_fixing(&ctx, 0, false, 1.0));
        // probing implications on a general integer or continuous column: not a clique row
        assert!(!s2_clique_change_is_fixing(&ctx, 1, true, 2.0));
        assert!(!s2_clique_change_is_fixing(&ctx, 1, true, 0.0));
        assert!(!s2_clique_change_is_fixing(&ctx, 1, false, 1.0));
        assert!(!s2_clique_change_is_fixing(&ctx, 2, true, 0.0));
    }

    #[test]
    fn s2_reason_clique_rows() {
        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };

        // ub change on y, val=1: -x - y >= -1
        clique_reason_row(0, 1, 1, true, &mut row);
        assert_eq!(row.vals, vec![-1.0, -1.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, -1.0);

        // ub change on y, val=0: x - y >= 0
        clique_reason_row(0, 0, 1, true, &mut row);
        assert_eq!(row.vals, vec![1.0, -1.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, 0.0);

        // lb change on y, val=1: -x + y >= 0
        clique_reason_row(0, 1, 1, false, &mut row);
        assert_eq!(row.vals, vec![-1.0, 1.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, 0.0);

        // lb change on y, val=0: x + y >= 1
        clique_reason_row(0, 0, 1, false, &mut row);
        assert_eq!(row.vals, vec![1.0, 1.0]);
        assert_eq!(row.inds, vec![0, 1]);
        assert_eq!(row.lhs, 1.0);
    }

    #[test]
    fn s2_reason_nogood() {
        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };

        let entries = vec![(0, false, 1.0), (1, true, 0.0)];
        let is_bin = |_c: HighsInt| true;
        assert!(nogood_reason_row(&entries, is_bin, &mut row));
        assert_eq!(row.vals, vec![-1.0, 1.0]);
        assert_eq!(row.lhs, 0.0);

        let is_bin_false = |col: HighsInt| col != 0;
        assert!(!nogood_reason_row(&entries, is_bin_false, &mut row));
    }

    #[test]
    fn s2_reason_proof() {
        let mut ms = tiny_mip();
        let mut dom = HighsDomain::new(&ms);
        let globaldom = HighsDomain::new(&ms);
        let mut ctx = S2DomainCtx {
            ms: &ms,
            localdom: &mut dom,
            globaldom: &globaldom,
        };

        let mut row = S2ConflictRow {
            vals: vec![],
            inds: vec![],
            lhs: 0.0,
            slack: 0.0,
            coefquotient: 0.0,
        };
        let inds = vec![0, 1];
        let vals = vec![1.0, 1.0];
        let upto = s2_initial_conflict_row(&mut ctx, Some((&inds, &vals, 2, 4.0)), &mut row);
        assert_eq!(upto, Some(kHighsIInf));
        assert_eq!(row.vals, vec![-1.0, -1.0]);
        assert_eq!(row.lhs, -4.0);
    }

    /// neos-960392: an objective-propagation reason
    /// row holds exact zeros. Kept in `inds`, they made linearCombRows append those columns a
    /// second time, and the resulting conflict cut had duplicate indices.
    #[test]
    fn s2_reason_zero_coefs_dropped() {
        use super::super::S2TestCtx;
        use super::super::row::linear_comb_rows;
        let ctx = S2TestCtx {
            glb: vec![0.0; 4],
            gub: vec![1.0, 1.0, 1.0, kHighsInf],
            integral: vec![true; 4],
            num_integral: 4,
            stack: vec![],
        };
        let reason = |vals: Vec<f64>| S2ConflictRow {
            inds: vec![0, 1, 2, 3],
            vals,
            lhs: -2.0,
            ..Default::default()
        };

        // Without the cleanup, linear_comb_rows duplicates the zero-coefficient columns.
        let mut raw = reason(vec![-1.0, 0.0, -1.0, 0.0]);
        let other = S2ConflictRow {
            inds: vec![1, 3],
            vals: vec![0.0, 1.0, 0.0, 1.0],
            lhs: 0.0,
            ..Default::default()
        };
        linear_comb_rows(&mut raw, &other, 1.0);
        assert!(!s2_inds_unique(&raw.inds));

        // Exact zeros are dropped and leave lhs unchanged.
        let mut row = reason(vec![-1.0, 0.0, -1.0, 0.0]);
        assert!(s2_drop_zero_coefs(&ctx, &mut row));
        let mut inds = row.inds.clone();
        inds.sort_unstable();
        assert_eq!(inds, vec![0, 2]);
        assert_eq!(row.lhs, -2.0);
        linear_comb_rows(&mut row, &other, 1.0);
        assert!(s2_inds_unique(&row.inds));
        assert_eq!(row.nnz(), 4);

        // A tiny coefficient is weakened out with its global bound: lhs -= a*ub.
        let mut row = reason(vec![-1.0, 1e-10, -1.0, 0.0]);
        assert!(s2_drop_zero_coefs(&ctx, &mut row));
        assert_eq!(row.nnz(), 2);
        assert_eq!(row.lhs, -2.0 - 1e-10);

        // A tiny coefficient on a column with an infinite bound rejects the row.
        let mut row = reason(vec![-1.0, 0.0, -1.0, 1e-10]);
        assert!(!s2_drop_zero_coefs(&ctx, &mut row));
    }
}
