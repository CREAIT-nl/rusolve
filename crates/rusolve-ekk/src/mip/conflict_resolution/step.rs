//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c (rescaleAndResolve 2188-2265, reduceReason 2270-2302,
//! slackReducingContinuousBdchgQueue 1730-1791, executeResolutionStep 2649-2812).

use std::collections::BinaryHeap;

use super::mir::mir_reduction;
use super::reason::get_reason_row;
use super::reduce::{complemented_mir_lhs, is_binary_conflict_row};
use super::row::{
    compute_scale_reason, get_quot_largest_smallest_coef, linear_comb_rows,
    weaken_continuous_vars_conflict_row,
};
use super::*;
use crate::prelude::*;

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::rescaleAndResolve (lines 2188-2265)
pub fn rescale_and_resolve(
    row1: &S2ConflictRow,
    row2: &S2ConflictRow,
    resolved: &mut S2ConflictRow,
    residx: HighsInt,
) -> bool {
    // scale = computeScaleReason(set, row1, row2, residx)
    let scale = compute_scale_reason(row1, row2, residx);

    // if( SCIPsetIsGE(set, scale, set->conf_maxcoefquot) )
    // {
    //    conflict->nreslargecoefs++;
    //    return SCIP_OKAY;
    // }
    // Note: conflict->nreslargecoefs++ is omitted; S2Stats has no nreslargecoefs field.
    if s2_is_ge(scale, S2_MAXCOEFQUOT) {
        return false;
    }

    // SCIP_CALL( conflictRowReplace(resolvedrow, blkmem, row1) );
    resolved.replace_from(row1);

    // newsize = resolvedrow->nnz + row2->nnz;
    // if( resolvedrow->size < newsize ) { ... }
    // (In Rust, replace_from already sized vals; linear_comb_rows handles inds growth.)

    // linearCombRows(set, resolvedrow, row2, scale);
    linear_comb_rows(resolved, row2, scale);

    // newnnz = resolvedrow->nnz;
    // largestcoef = -SCIPsetInfinity(set);
    // smallestcoef = SCIPsetInfinity(set);
    // for( i = newnnz - 1; i >= 0 ; i-- ) { ... }
    // (get_quot_largest_smallest_coef computes the same quotient in one pass.)

    // resolvedrow->coefquotient = (resolvedrow->nnz > 0) ? fabs(largestcoef / smallestcoef) : 0.0;
    resolved.coefquotient =
        get_quot_largest_smallest_coef(&resolved.inds, &resolved.vals, resolved.nnz());

    // if( SCIPsetIsGT(set, resolvedrow->coefquotient, set->conf_maxcoefquot) )
    // {
    //    conflict->nreslargecoefs++;
    //    return SCIP_OKAY;
    // }
    // Note: conflict->nreslargecoefs++ is omitted; S2Stats has no nreslargecoefs field.
    if s2_is_gt(resolved.coefquotient, S2_MAXCOEFQUOT) {
        return false;
    }

    // *success = TRUE;
    true
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::reduceReason (lines 2270-2302)
pub fn reduce_reason<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    fixbounds: &[f64],
    fixsides: &[i8],
    pos: HighsInt,
    residx: HighsInt,
) {
    // coefinrow = fabs(rowtoreduce->vals[residx]);
    let coefinrow = row.vals[residx as usize].abs();

    // if( isBinaryConflictRow(set, vars, rowtoreduce) )
    //    SCIP_CALL( ComplementedMirLhs(set, vars, rowtoreduce, fixsides, currbdchgidx, residx, coefinrow) );
    // else
    //    SCIP_CALL( MirReduction(set, blkmem, vars, nvars, rowtoreduce, currbdchginfo, residx, coefinrow) );
    if is_binary_conflict_row(ctx, row) {
        complemented_mir_lhs(ctx, row, fixsides, pos, residx, coefinrow);
    } else {
        let _ = mir_reduction(ctx, row, pos, residx, coefinrow);
    }

    // SCIP_CALL( computeSlack(set, vars, rowtoreduce, currbdchginfo, fixbounds, fixsides) );
    compute_slack(ctx, row, pos, Some((fixbounds, fixsides)));

    // assert(!SCIPsetIsZero(set, rowtoreduce->vals[residx]));
    debug_assert!(!s2_is_zero(row.vals[residx as usize]));
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::slackReducingContinuousBdchgQueue (lines 1730-1791)
pub fn slack_reducing_continuous_queue<C: S2Context>(
    ctx: &C,
    q: &mut BinaryHeap<HighsInt>,
    row: &S2ConflictRow,
    inferpos: HighsInt,
) {
    // make sure that the continuous bound change queue is empty
    q.clear();

    // for( i = 0; i < row->nnz; i++ )
    for i in 0..row.nnz() {
        // v = row->inds[i];
        let v = row.inds[i];

        // if( SCIPvarIsIntegral(vars[v]) ) continue;
        if ctx.is_integral(v) {
            continue;
        }

        // coef = row->vals[v];
        let coef = row.vals[v as usize];

        // if( coef > 0.0 )
        if coef > 0.0 {
            // SCIP_BDCHGINFO* bdchginfo = SCIPvarGetUbchgInfo(vars[v], inferbdchgidx, FALSE);
            // after = FALSE means upto = inferpos - 1
            let (_bound, p) = ctx.ub_upto(v, inferpos - 1);
            // if( bdchginfo != NULL ) { SCIP_CALL( SCIPpqueueInsert(...) ); }
            if p != -1 {
                q.push(p);
            }
        } else {
            // SCIP_BDCHGINFO* bdchginfo = SCIPvarGetLbchgInfo(vars[v], inferbdchgidx, FALSE);
            let (_bound, p) = ctx.lb_upto(v, inferpos - 1);
            if p != -1 {
                q.push(p);
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::executeResolutionStep (lines 2649-2812)
pub fn execute_resolution_step<C: S2Context>(
    ctx: &mut C,
    ws: &mut HighsConflictResolution,
    pos: HighsInt,
    residx: HighsInt,
    maxsize: HighsInt,
) -> bool {
    // conflictrow = conflict->conflictrow;
    // reasonrow = conflict->reasonrow;
    // reducedreasonrow = conflict->reducedreasonrow;
    // resolvedconflictrow = conflict->resolvedconflictrow;
    //
    // In Rust we use std::mem::take to move rows out of ws for mutable calls,
    // then put them back, keeping the C call order.

    let mut conflictrow = std::mem::take(&mut ws.conflictrow);
    let mut reasonrow = std::mem::take(&mut ws.reasonrow);
    let mut reducedreasonrow = std::mem::take(&mut ws.reducedreasonrow);
    let mut resolvedconflictrow = std::mem::take(&mut ws.resolvedconflictrow);

    // nvars = SCIPgetNVars(set->scip);
    let nvars = ctx.num_col();

    // assert(!SCIPsetIsZero(set, reasonrow->vals[residx]));
    debug_assert!(!s2_is_zero(reasonrow.vals[residx as usize]));

    // fixbounds and fixsides are borrowed from ws; we take them out too.
    let mut fixbounds = std::mem::take(&mut ws.fixbounds);
    let mut fixsides = std::mem::take(&mut ws.fixsides);

    // try to resolve without reducing the reason row
    // SCIP_CALL( rescaleAndResolve(set, conflict, conflictrow, reasonrow, resolvedconflictrow, blkmem,
    //                      residx, successresolution) );
    let mut success =
        rescale_and_resolve(&conflictrow, &reasonrow, &mut resolvedconflictrow, residx);

    // if( !(*successresolution) ) return SCIP_OKAY;
    if !success {
        ws.conflictrow = conflictrow;
        ws.reasonrow = reasonrow;
        ws.reducedreasonrow = reducedreasonrow;
        ws.resolvedconflictrow = resolvedconflictrow;
        ws.fixbounds = fixbounds;
        ws.fixsides = fixsides;
        return success;
    }

    // SCIP_CALL( computeSlack(set, vars, resolvedconflictrow, currbdchginfo, fixbounds, fixsides) );
    compute_slack(
        ctx,
        &mut resolvedconflictrow,
        pos,
        Some((&fixbounds, &fixsides)),
    );

    // the C line `if( set->conf_reduction == 'o' ) return` is OMITTED:
    // conf/reduction defaults to 'm' (set.c:184) and no constant for it exists in D7.
    s2_trace_row("step.rescale1", &resolvedconflictrow, &ws.debug_sol);

    // if( SCIPsetIsGE(set, resolvedconflictrow->slack, 0.0) )
    if s2_is_ge(resolvedconflictrow.slack, 0.0) {
        // copy the original reason here
        // SCIP_CALL( conflictRowReplace(reducedreasonrow, blkmem, reasonrow) );
        reducedreasonrow.replace_from(&reasonrow);
        s2_trace_row("step.reduce.copy", &reducedreasonrow, &ws.debug_sol);

        // apply reduction to the reason row
        // SCIP_CALL( reduceReason(set, blkmem, vars, nvars, reducedreasonrow, currbdchginfo, residx, fixbounds, fixsides) );
        reduce_reason(
            ctx,
            &mut reducedreasonrow,
            &fixbounds,
            &fixsides,
            pos,
            residx,
        );
        s2_trace_row("step.reduce.mir", &reducedreasonrow, &ws.debug_sol);

        // after reduction resolve again
        // SCIP_CALL( rescaleAndResolve(set, conflict, conflictrow, reducedreasonrow, resolvedconflictrow, blkmem,
        //                         residx, successresolution) );
        success = rescale_and_resolve(
            &conflictrow,
            &reducedreasonrow,
            &mut resolvedconflictrow,
            residx,
        );

        // if( !(*successresolution) ) return SCIP_OKAY;
        if !success {
            ws.conflictrow = conflictrow;
            ws.reasonrow = reasonrow;
            ws.reducedreasonrow = reducedreasonrow;
            ws.resolvedconflictrow = resolvedconflictrow;
            ws.fixbounds = fixbounds;
            ws.fixsides = fixsides;
            return success;
        }

        // SCIP_CALL( computeSlack(set, vars, resolvedconflictrow, currbdchginfo, fixbounds, fixsides) );
        compute_slack(
            ctx,
            &mut resolvedconflictrow,
            pos,
            Some((&fixbounds, &fixsides)),
        );
        s2_trace_row("step.rescale2", &resolvedconflictrow, &ws.debug_sol);

        // mixed binary case reduction
        // if( set->conf_mbreduction && SCIPsetIsGE(set, resolvedconflictrow->slack, 0.0) )
        if S2_MBREDUCTION && s2_is_ge(resolvedconflictrow.slack, 0.0) {
            // now both the reducedreasonrow and reasonrow are the rows before the reduction
            // SCIP_CALL( conflictRowReplace(reducedreasonrow, blkmem, reasonrow) );
            reducedreasonrow.replace_from(&reasonrow);
            s2_trace_row("step.mb.copy", &reducedreasonrow, &ws.debug_sol);

            // coefresidx = reducedreasonrow->vals[residx];
            let coefresidx = reducedreasonrow.vals[residx as usize];

            // In the following the reducedreasonrow is the aggregation of the different reasonrows of continuous bound changes
            // SCIP_CALL( slackReducingContinuousBdchgQueue(conflict, vars, reducedreasonrow, SCIPbdchginfoGetIdx(currbdchginfo)) );
            slack_reducing_continuous_queue(ctx, &mut ws.contqueue, &reducedreasonrow, pos);

            // if( SCIPpqueueNElems(conflict->continuousbdchgqueue) > 0 ) { ... } else return SCIP_OKAY;
            if ws.contqueue.is_empty() {
                ws.conflictrow = conflictrow;
                ws.reasonrow = reasonrow;
                ws.reducedreasonrow = reducedreasonrow;
                ws.resolvedconflictrow = resolvedconflictrow;
                ws.fixbounds = fixbounds;
                ws.fixsides = fixsides;
                return success;
            }

            // int nresconts = 0;
            let mut nresconts: HighsInt = 0;

            // while( SCIPpqueueNElems(conflict->continuousbdchgqueue) > 0 )
            while !ws.contqueue.is_empty() {
                // nresconts++;
                nresconts += 1;

                // if( nresconts > maxsize / 2 || reducedreasonrow->nnz > maxsize )
                if nresconts > maxsize / 2 || reducedreasonrow.nnz() as HighsInt > maxsize {
                    // (*successresolution) = FALSE;
                    success = false;
                    // conflict->nreslongconfs++;
                    ws.stats.nreslongconfs += 1;
                    ws.conflictrow = conflictrow;
                    ws.reasonrow = reasonrow;
                    ws.reducedreasonrow = reducedreasonrow;
                    ws.resolvedconflictrow = resolvedconflictrow;
                    ws.fixbounds = fixbounds;
                    ws.fixsides = fixsides;
                    return success;
                }

                // continuousbdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->continuousbdchgqueue));
                let contpos = ws.contqueue.pop().unwrap();

                // int varidx = SCIPvarGetProbindex(SCIPbdchginfoGetVar(continuousbdchginfo));
                let varidx = ctx.change(contpos).col;

                // conflictRowClear(blkmem, reasonrow, nvars);
                reasonrow.clear(nvars);

                // SCIP_CALL( getReasonRow(blkmem, vars, set, continuousbdchginfo, reasonrow, varidx, fixbounds, fixsides,
                //                         &successgetreason) );
                let successgetreason =
                    get_reason_row(ctx, contpos, &mut reasonrow, varidx, &fixbounds, &fixsides);
                s2_trace_row("step.cont.reason", &reasonrow, &ws.debug_sol);

                // if( !successgetreason ) { (*successresolution) = FALSE; return SCIP_OKAY; }
                if !successgetreason {
                    success = false;
                    ws.conflictrow = conflictrow;
                    ws.reasonrow = reasonrow;
                    ws.reducedreasonrow = reducedreasonrow;
                    ws.resolvedconflictrow = resolvedconflictrow;
                    ws.fixbounds = fixbounds;
                    ws.fixsides = fixsides;
                    return success;
                }

                // resolve the continuous variable
                // SCIP_CALL( rescaleAndResolve(set, conflict, reducedreasonrow, reasonrow, resolvedconflictrow, blkmem,
                //                      varidx, successresolution) );
                success = rescale_and_resolve(
                    &reducedreasonrow,
                    &reasonrow,
                    &mut resolvedconflictrow,
                    varidx,
                );

                // if( !(*successresolution) ) return SCIP_OKAY;
                if !success {
                    ws.conflictrow = conflictrow;
                    ws.reasonrow = reasonrow;
                    ws.reducedreasonrow = reducedreasonrow;
                    ws.resolvedconflictrow = resolvedconflictrow;
                    ws.fixbounds = fixbounds;
                    ws.fixsides = fixsides;
                    return success;
                }

                // if the coefficient of the variable we are resolving changed sign, then we have a new conflict constraint
                // if( SCIPsetIsLE(set, coefresidx * resolvedconflictrow->vals[residx], 0.0) )
                if s2_is_le(coefresidx * resolvedconflictrow.vals[residx as usize], 0.0) {
                    // #ifndef SCIP_DEBUG
                    // SCIP_CALL( computeSlack(set, vars, resolvedconflictrow, currbdchginfo, fixbounds, fixsides) );
                    // assert(SCIPsetIsLT(set, resolvedconflictrow->slack, 0.0) || !isBinaryConflictRow(set, vars, reducedreasonrow));
                    // #endif
                    compute_slack(
                        ctx,
                        &mut resolvedconflictrow,
                        pos,
                        Some((&fixbounds, &fixsides)),
                    );
                    debug_assert!(
                        s2_is_lt(resolvedconflictrow.slack, 0.0)
                            || !is_binary_conflict_row(ctx, &reducedreasonrow)
                    );
                    ws.conflictrow = conflictrow;
                    ws.reasonrow = reasonrow;
                    ws.reducedreasonrow = reducedreasonrow;
                    ws.resolvedconflictrow = resolvedconflictrow;
                    ws.fixbounds = fixbounds;
                    ws.fixsides = fixsides;
                    return success;
                }

                // SCIP_CALL( conflictRowReplace(reducedreasonrow, blkmem, resolvedconflictrow) );
                reducedreasonrow.replace_from(&resolvedconflictrow);
                s2_trace_row("step.cont.resolved", &reducedreasonrow, &ws.debug_sol);

                // SCIP_CALL( slackReducingContinuousBdchgQueue(conflict, vars, reducedreasonrow, SCIPbdchginfoGetIdx(continuousbdchginfo)) );
                slack_reducing_continuous_queue(ctx, &mut ws.contqueue, &reducedreasonrow, contpos);
            }

            // after the loop we can weaken all continuous variables with their global bounds
            // SCIP_CALL( weakenContinuousVarsConflictRow(reducedreasonrow, set, vars, residx) );
            weaken_continuous_vars_conflict_row(ctx, &mut reducedreasonrow, residx);
            s2_trace_row("step.mb.weakencont", &reducedreasonrow, &ws.debug_sol);

            // after resolving the continuous variables we resolve once more if the constraint is binary
            // if( nresconts > 0 )
            if nresconts > 0 {
                // apply reduction to the reason row
                // SCIP_CALL( reduceReason(set, blkmem, vars, nvars, reducedreasonrow, currbdchginfo, residx, fixbounds, fixsides) );
                reduce_reason(
                    ctx,
                    &mut reducedreasonrow,
                    &fixbounds,
                    &fixsides,
                    pos,
                    residx,
                );
                s2_trace_row("step.mb.reduce", &reducedreasonrow, &ws.debug_sol);

                // after all reductions resolve one final time
                // SCIP_CALL( rescaleAndResolve(set, conflict, conflictrow, reducedreasonrow, resolvedconflictrow, blkmem,
                //                      residx, successresolution) );
                success = rescale_and_resolve(
                    &conflictrow,
                    &reducedreasonrow,
                    &mut resolvedconflictrow,
                    residx,
                );

                // SCIP_CALL( computeSlack(set, vars, resolvedconflictrow, currbdchginfo, fixbounds, fixsides) );
                compute_slack(
                    ctx,
                    &mut resolvedconflictrow,
                    pos,
                    Some((&fixbounds, &fixsides)),
                );
            }
        }
    }

    // put rows back
    ws.conflictrow = conflictrow;
    ws.reasonrow = reasonrow;
    ws.reducedreasonrow = reducedreasonrow;
    ws.resolvedconflictrow = resolvedconflictrow;
    ws.fixbounds = fixbounds;
    ws.fixsides = fixsides;

    success
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn s2_step_rescale_basic() {
        // conflict row: x1 + x2 >= 1, inds [1,2], vals [0,1,1], lhs 1.0
        let mut conflict = S2ConflictRow::default();
        conflict.ensure_dim(3);
        conflict.inds = vec![1, 2];
        conflict.vals[1] = 1.0;
        conflict.vals[2] = 1.0;
        conflict.lhs = 1.0;

        // reason row: -x0 - x2 >= -1, inds [0,2], vals [-1,0,-1], lhs -1.0
        let mut reason = S2ConflictRow::default();
        reason.ensure_dim(3);
        reason.inds = vec![0, 2];
        reason.vals[0] = -1.0;
        reason.vals[2] = -1.0;
        reason.lhs = -1.0;

        let mut resolved = S2ConflictRow::default();
        resolved.ensure_dim(3);

        let result = rescale_and_resolve(&conflict, &reason, &mut resolved, 2);
        assert!(result);

        // scale = |coefconf/coefreas| = |1/-1| = 1
        // resolved = conflict + 1 * reason:
        //   x1 stays 1, x2 cancels to 0 and is removed, x0 becomes -1
        //   lhs = 1 + (-1) = 0
        // Inds order from replace_from (inds [1,2]) then linear_comb_rows pushes x0:
        //   replace_from copies inds [1,2], then linear_comb_rows adds x0 (inds.push(0))
        //   x2 cancels and is removed, leaving inds [1, 0]
        assert_eq!(resolved.inds, vec![1, 0]);
        assert_eq!(resolved.vals[1], 1.0);
        assert_eq!(resolved.vals[0], -1.0);
        assert_eq!(resolved.lhs, 0.0);
        assert_eq!(resolved.coefquotient, 1.0);
    }

    #[test]
    fn s2_step_rescale_limit() {
        // conflict: x0 with 1e7, x2 with 1, inds [0,2], lhs 1.0
        let mut conflict = S2ConflictRow::default();
        conflict.ensure_dim(3);
        conflict.inds = vec![0, 2];
        conflict.vals[0] = 1e7;
        conflict.vals[2] = 1.0;
        conflict.lhs = 1.0;

        // reason: x0 with -1, inds [0], lhs -1.0
        let mut reason = S2ConflictRow::default();
        reason.ensure_dim(3);
        reason.inds = vec![0];
        reason.vals[0] = -1.0;
        reason.lhs = -1.0;

        let mut resolved = S2ConflictRow::default();
        resolved.ensure_dim(3);

        // scale = |1e7 / -1| = 1e7 >= S2_MAXCOEFQUOT (1e6)
        let result = rescale_and_resolve(&conflict, &reason, &mut resolved, 0);
        assert!(!result);
    }

    #[test]
    fn s2_step_reduce_binary() {
        // S2TestCtx: 6 binary columns, num_integral 6
        let ctx = S2TestCtx {
            glb: vec![0.0; 6],
            gub: vec![1.0; 6],
            integral: vec![true; 6],
            num_integral: 6,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: true,
                    reason: None,
                },
                S2TestChg {
                    col: 2,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: true,
                    reason: None,
                },
                S2TestChg {
                    col: 4,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: true,
                    reason: None,
                },
                S2TestChg {
                    col: 5,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: true,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 1.0,
                    depth: 1,
                    branching: false,
                    reason: Some((vec![(0, 2.0), (1, 1.0), (2, 1.0)], 1.0)),
                },
            ],
        };

        // Conflict C: -2x0 + x4 + x5 >= -1, inds [0,4,5]
        let mut conflict = S2ConflictRow::default();
        conflict.ensure_dim(6);
        conflict.inds = vec![0, 4, 5];
        conflict.vals[0] = -2.0;
        conflict.vals[4] = 1.0;
        conflict.vals[5] = 1.0;
        conflict.lhs = -1.0;

        // Reason R: 2x0 + x1 + x2 >= 1, inds [0,1,2]
        let mut reason = S2ConflictRow::default();
        reason.ensure_dim(6);
        reason.inds = vec![0, 1, 2];
        reason.vals[0] = 2.0;
        reason.vals[1] = 1.0;
        reason.vals[2] = 1.0;
        reason.lhs = 1.0;

        let mut ws = HighsConflictResolution::default();
        ws.ensure_dim(6);
        ws.conflictrow = conflict;
        ws.reasonrow = reason;

        let mut ctx = ctx;
        let result = execute_resolution_step(&mut ctx, &mut ws, 4, 0, 10);

        assert!(result);

        // Trace: the first resolvent C + R has slack 0 (all used columns at ub 0), so the
        // reason is reduced by complemented_mir_lhs (divisor 2) to x0 + x1 + x2 >= 1.
        // The second resolve uses scale 2. linear_comb_rows then swap_removes x0 (zero), giving
        // inds [2,4,5,1] with x1 = 2, x2 = 2, x4 = 1, x5 = 1, lhs 1.0, slack -1.0.
        // The mixed-binary branch is not entered because slack < 0.
        assert_eq!(ws.resolvedconflictrow.vals[1], 2.0);
        assert_eq!(ws.resolvedconflictrow.vals[2], 2.0);
        assert_eq!(ws.resolvedconflictrow.vals[4], 1.0);
        assert_eq!(ws.resolvedconflictrow.vals[5], 1.0);
        assert_eq!(ws.resolvedconflictrow.vals[0], 0.0);
        assert_eq!(ws.resolvedconflictrow.lhs, 1.0);
        assert_eq!(ws.resolvedconflictrow.slack, -1.0);
        assert_eq!(ws.resolvedconflictrow.coefquotient, 2.0);
        assert_eq!(ws.resolvedconflictrow.inds, vec![2, 4, 5, 1]);
    }

    #[test]
    fn s2_step_cont_queue() {
        // Row: x0 + y >= 1, inds [0,1], vals [1,1], lhs 1.0
        let mut row = S2ConflictRow::default();
        row.ensure_dim(2);
        row.inds = vec![0, 1];
        row.vals[0] = 1.0;
        row.vals[1] = 1.0;
        row.lhs = 1.0;

        // Columns: col 0 integral, col 1 (y) continuous
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![1.0, 5.0],
            integral: vec![true, false],
            num_integral: 1,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 2.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 1.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
            ],
        };

        // inferpos 3: ub_upto(col 1, 2) -> latest ub change at or before pos 2 is pos 2 (y ub 0).
        //   So p = 2. q = [2].
        let mut q = BinaryHeap::new();
        slack_reducing_continuous_queue(&ctx, &mut q, &row, 3);
        assert_eq!(q.into_sorted_vec(), vec![2]);

        // inferpos 2: ub_upto(col 1, 1) -> latest ub change at or before pos 1 is pos 0 (y ub 2).
        //   So p = 0. q = [0].
        let mut q2 = BinaryHeap::new();
        slack_reducing_continuous_queue(&ctx, &mut q2, &row, 2);
        assert_eq!(q2.into_sorted_vec(), vec![0]);
    }
}
