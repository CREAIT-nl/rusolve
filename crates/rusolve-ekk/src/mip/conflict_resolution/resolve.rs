//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictCalcResMaxsize (lines 487-500),
//! conflictInsertConflictRow (208-230), getConflictRow (2572-2645), markBdchgAsFixed (2819-2904),
//! conflictAnalyzeResolution (2960-3333; addConflictRows call excluded).

use super::queue::update_bdchg_queue;
use super::reason::get_reason_row;
use super::row::{get_quot_largest_smallest_coef, tighten_coefs_traced, weaken_conflict_row};
use super::step::execute_resolution_step;
use super::*;
use crate::prelude::*;
use crate::util::highs_utils::cpp_max;

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictCalcResMaxsize (lines 487-500)
pub fn conflict_calc_res_maxsize<C: S2Context>(ctx: &C) -> HighsInt {
    cpp_max(
        (S2_MAXVARSFRACRES * ctx.num_integral() as f64) as HighsInt,
        S2_MINMAXVARS,
    )
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictInsertConflictRow (lines 208-230)
fn conflict_insert_conflict_row(ws: &mut HighsConflictResolution, row: &S2ConflictRow) {
    ws.conflictrows.push(row.clone());
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::getConflictRow (lines 2572-2645)
fn get_conflict_row<C: S2Context>(
    ctx: &C,
    ws: &mut HighsConflictResolution,
    initial: &S2ConflictRow,
    pos: HighsInt,
    maxsize: HighsInt,
) -> bool {
    // C conflictRowFromLpRow is replaced by this copy; the initial row is always given here.
    ws.conflictrow.replace_from(initial);
    s2_trace_row("get_conflict_row.copy", &ws.conflictrow, &ws.debug_sol);

    // if the conflict row is too large, we try to weaken it
    if ws.conflictrow.nnz() as HighsInt > maxsize {
        weaken_conflict_row(ctx, &mut ws.conflictrow, pos, None);
        s2_trace_row("get_conflict_row.weaken", &ws.conflictrow, &ws.debug_sol);

        if ws.conflictrow.nnz() as HighsInt > maxsize {
            ws.stats.nreslongconfs += 1;
            return false;
        }
    }

    // set the slack
    compute_slack(ctx, &mut ws.conflictrow, pos, None);

    // The conflict row might not be infeasible.
    if s2_is_ge(ws.conflictrow.slack, 0.0) {
        return false;
    }

    // check once more if the initial conflict is too long
    if ws.conflictrow.nnz() as HighsInt > maxsize {
        ws.stats.nreslongconfs += 1;
        return false;
    }

    debug_assert!(s2_is_lt(ws.conflictrow.slack, 0.0));
    true
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::markBdchgAsFixed (lines 2819-2904)
fn mark_bdchg_as_fixed<C: S2Context>(
    ctx: &C,
    ws: &mut HighsConflictResolution,
    pos: &mut HighsInt,
    depth: &mut HighsInt,
    nressteps: HighsInt,
) -> bool {
    let chg = ctx.change(*pos);

    // The C test `if SCIPbdchginfoGetChgtype == CONSINFER && !SCIPconsIsGlobal` is always
    // false here: all HiGHS reasons are global.

    if !ws.resqueue.exists_resolvable(ctx) {
        return false;
    }

    let varidx = chg.col as usize;

    // If a bound for this variable has already been ignored, abort
    if ws.fixsides[varidx] != 0 {
        return false;
    }

    // Record the ignored bound change: 1 for upper, -1 for lower
    ws.fixsides[varidx] = if chg.is_upper { 1 } else { -1 };
    ws.fixbounds[varidx] = chg.newbound;

    // Reset conflict bounds to allow weaker bounds to be added later
    if !chg.is_upper {
        ws.resqueue.conflictvarslbs[varidx] = f64::MIN;
    } else {
        ws.resqueue.conflictvarsubs[varidx] = f64::MAX;
    }

    if nressteps == 0 {
        // C: updateBdchgQueue with after=FALSE, so upto = pos - 1
        update_bdchg_queue(ctx, &mut ws.resqueue, &ws.conflictrow, *pos - 1);
    }

    let first = ws.resqueue.first_cand(ctx, &ws.conflictrow, false);
    debug_assert!(first.is_some());

    let popped = match ws.resqueue.remove_cand() {
        Some(p) => p,
        None => {
            debug_assert!(false, "remove_cand returned None after first_cand was Some");
            return false;
        }
    };

    *pos = popped;

    // If no resolution step has been applied and the bound change is a branching decision,
    // update the last bound change depth.
    if nressteps == 0 && chg.branching {
        *depth = ctx.change(*pos).depth;
    } else if !S2_FIXANDCONTINUE {
        return false;
    }

    // #ifndef SCIP_DEBUG: computeSlack is ported
    compute_slack(
        ctx,
        &mut ws.conflictrow,
        *pos,
        Some((&ws.fixbounds, &ws.fixsides)),
    );
    debug_assert!(s2_is_lt(ws.conflictrow.slack, 0.0));

    true
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictAnalyzeResolution (lines 2960-3333; addConflictRows call excluded)
fn run_resolution<C: S2Context>(
    ctx: &mut C,
    ws: &mut HighsConflictResolution,
    initial: &S2ConflictRow,
    upto: HighsInt,
) {
    let nvars = ctx.num_col();
    s2_trace_row("initial", initial, &ws.debug_sol);
    let mut lastuipdepth: HighsInt = -1;

    // clear the conflict, reason, resolved conflict rows
    ws.conflictrow.clear(nvars);
    ws.resolvedconflictrow.clear(nvars);
    ws.reasonrow.clear(nvars);
    ws.reducedreasonrow.clear(nvars);

    // reset the queues: SCIP clears both at the end of the previous call
    // (SCIPconflictAnalyzeResolution 3400-3401) and resets the bounds in
    // conflict_graphanalysis.c:3435-3451; doing it here is equivalent.
    ws.resqueue.reset(nvars);
    ws.contqueue.clear();

    // seed the queue (D4/D5 step 1)
    update_bdchg_queue(ctx, &mut ws.resqueue, initial, upto);

    // last bound change that led to infeasibility
    let Some(_first) = ws.resqueue.first_cand(ctx, &ws.conflictrow, true) else {
        return;
    };

    // if no bound change exists or none of them was inferred by a resolvable
    // constraint then we terminate
    if !ws.resqueue.exists_resolvable(ctx) {
        return;
    }

    // remove the last bound change
    let Some(mut pos) = ws.resqueue.remove_cand() else {
        return;
    };

    let chg = ctx.change(pos);
    let mut depth = chg.depth;

    // calculate the maximal size of each accepted conflict set
    let maxsize: HighsInt = 2 * conflict_calc_res_maxsize(ctx);

    // sets the initial conflictrow
    if !get_conflict_row(ctx, ws, initial, pos, maxsize) {
        return;
    }

    // before we resolve any bound change the resolved conflict row is the same as the conflict row
    ws.resolvedconflictrow.replace_from(&ws.conflictrow);

    // Apply coefficient tightening to the conflict constraint should never hurt
    let mut nnz = ws.conflictrow.nnz();
    let nchgcoefs = tighten_coefs_traced(
        ctx,
        &mut ws.conflictrow.vals,
        &mut ws.conflictrow.inds,
        &mut nnz,
        &mut ws.conflictrow.lhs,
        &ws.debug_sol,
    );
    ws.conflictrow.inds.truncate(nnz);
    s2_trace_row("setup.tighten", &ws.conflictrow, &ws.debug_sol);

    if nchgcoefs > 0 {
        // remove variables with zero coefficients from the conflict row
        remove_zero_vars(&mut ws.conflictrow);
        compute_slack(ctx, &mut ws.conflictrow, pos, None);
        // the slack after coefficient tightening must also be negative
        debug_assert!(s2_is_lt(ws.conflictrow.slack, 0.0));
    }

    ws.conflictrow.coefquotient = get_quot_largest_smallest_coef(
        &ws.conflictrow.inds,
        &ws.conflictrow.vals,
        ws.conflictrow.nnz(),
    );

    // initialize indices and bounds for the unresolvable bound changes
    ws.fixsides.clear();
    ws.fixsides.resize(nvars, 0);
    ws.fixbounds.clear();
    ws.fixbounds.resize(nvars, 0.0);

    // main loop
    loop {
        let chg = ctx.change(pos);

        if !chg.resolvable {
            let nres = ws.nressteps;
            if !mark_bdchg_as_fixed(ctx, ws, &mut pos, &mut depth, nres) {
                break;
            }
            s2_trace_row("mark_fixed", &ws.conflictrow, &ws.debug_sol);
            // C re-reads bdchgidx and vartoresolve from the new bdchginfo;
            // in Rust these are computed fresh from ctx.change(pos) each iteration.
        } else {
            let residx = chg.col;

            ws.reasonrow.clear(nvars);

            // get reason row of the latest bdchginfo
            if !get_reason_row(
                ctx,
                pos,
                &mut ws.reasonrow,
                residx,
                &ws.fixbounds,
                &ws.fixsides,
            ) {
                break;
            }
            s2_trace_row("reasonrow", &ws.reasonrow, &ws.debug_sol);

            // call resolution
            let success = execute_resolution_step(ctx, ws, pos, residx, maxsize);

            if !success {
                break;
            }

            ws.conflictrow.replace_from(&ws.resolvedconflictrow);
            s2_trace_row("resolved", &ws.conflictrow, &ws.debug_sol);

            // we must reset the conflict lower and upper bound to be able to add weaker bounds later
            if !chg.is_upper {
                ws.resqueue.conflictvarslbs[residx as usize] = f64::MIN;
            } else {
                ws.resqueue.conflictvarsubs[residx as usize] = f64::MAX;
            }

            if ws.conflictrow.nnz() as HighsInt > maxsize {
                weaken_conflict_row(ctx, &mut ws.conflictrow, pos, Some(ws.fixsides.as_slice()));
                s2_trace_row("loop.weaken", &ws.conflictrow, &ws.debug_sol);
                if ws.conflictrow.nnz() as HighsInt > maxsize {
                    ws.stats.nreslongconfs += 1;
                    break;
                }
            }
            // C debug-only computeSlack after weakening is NOT ported

            // A positive slack means that the conflict row is not infeasible anymore.
            if s2_is_ge(ws.conflictrow.slack, 0.0) {
                break;
            }

            ws.nressteps += 1;

            // apply coefficient tightening to the resolved constraint should never hurt
            let mut nnz = ws.conflictrow.nnz();
            let nchgcoefs = tighten_coefs_traced(
                ctx,
                &mut ws.conflictrow.vals,
                &mut ws.conflictrow.inds,
                &mut nnz,
                &mut ws.conflictrow.lhs,
                &ws.debug_sol,
            );
            ws.conflictrow.inds.truncate(nnz);
            s2_trace_row("loop.tighten", &ws.conflictrow, &ws.debug_sol);

            if nchgcoefs > 0 {
                // remove variables with zero coefficients from the conflict row
                remove_zero_vars(&mut ws.conflictrow);
                // The new slack should always be less or equal to the old slack
                compute_slack(
                    ctx,
                    &mut ws.conflictrow,
                    pos,
                    Some((&ws.fixbounds, &ws.fixsides)),
                );
            }

            // if we reached this point the conflict constraint must have negative slack
            debug_assert!(s2_is_lt(ws.conflictrow.slack, 0.0));

            // C FALSE: after=false, so upto = pos - 1
            update_bdchg_queue(ctx, &mut ws.resqueue, &ws.conflictrow, pos - 1);

            // get the next bound change
            if ws
                .resqueue
                .first_cand(ctx, &ws.conflictrow, false)
                .is_none()
            {
                // reached UIP; add the previous conflict in the list of conflict rows
                let row = ws.conflictrow.clone();
                conflict_insert_conflict_row(ws, &row);
                break;
            }

            let Some(popped) = ws.resqueue.remove_cand() else {
                break;
            };
            pos = popped;

            let chg2 = ctx.change(pos);
            depth = chg2.depth;

            debug_assert!(!s2_is_zero(ws.conflictrow.vals[chg2.col as usize]));

            // get the bound change before bdchginfo
            let nextbdchg = ws.resqueue.first_cand(ctx, &ws.conflictrow, false);

            debug_assert!(nextbdchg.is_none_or(|p| ctx.change(p).depth <= depth));

            // when at an UIP add the previous conflict in the list of conflict rows
            if nextbdchg.is_none_or(|p| ctx.change(p).depth != depth) && chg2.depth != lastuipdepth
            {
                let row = ws.conflictrow.clone();
                conflict_insert_conflict_row(ws, &row);
                lastuipdepth = chg2.depth;
                ws.nfuips += 1;
                // stop after conf_resfuiplevels UIPs
                if S2_RESFUIPLEVELS > 0 && ws.nfuips >= S2_RESFUIPLEVELS {
                    break;
                }
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::conflictAnalyzeResolution (lines 2960-3333)
pub fn analyze_resolution<C: S2Context>(
    ctx: &mut C,
    ws: &mut HighsConflictResolution,
    initial: &S2ConflictRow,
    upto: HighsInt,
) -> usize {
    ws.conflictrows.clear();
    ws.nressteps = 0;
    ws.nfuips = 0;

    run_resolution(ctx, ws, initial, upto);

    ws.stats.nressteps += ws.nressteps as i64;
    ws.stats.rows_found += ws.conflictrows.len() as i64;

    ws.conflictrows.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::conflict_resolution::{S2TestChg, S2TestCtx};
    use crate::prelude::*;

    fn s2_res_ctx(n: usize, nint: usize, stack: Vec<S2TestChg>) -> S2TestCtx {
        let integral: Vec<bool> = (0..n).map(|i| i < nint).collect();
        let glb = vec![0.0; n];
        let gub = vec![1.0; n];
        S2TestCtx {
            glb,
            gub,
            integral,
            num_integral: nint,
            stack,
        }
    }

    fn s2_res_initial(n: usize, terms: &[(HighsInt, f64)], lhs: f64) -> S2ConflictRow {
        let mut row = S2ConflictRow::default();
        row.ensure_dim(n);
        row.inds = terms.iter().map(|&(col, _)| col).collect();
        for &(col, v) in terms {
            row.vals[col as usize] = v;
        }
        row.lhs = lhs;
        row
    }

    fn chain_stack() -> Vec<S2TestChg> {
        vec![
            S2TestChg {
                col: 0,
                is_upper: false,
                newbound: 1.0,
                depth: 1,
                branching: true,
                reason: None,
            },
            S2TestChg {
                col: 1,
                is_upper: true,
                newbound: 0.0,
                depth: 1,
                branching: false,
                reason: Some((vec![(0, -1.0), (1, -1.0)], -1.0)),
            },
            S2TestChg {
                col: 2,
                is_upper: true,
                newbound: 0.0,
                depth: 1,
                branching: false,
                reason: Some((vec![(0, -1.0), (2, -1.0)], -1.0)),
            },
        ]
    }

    #[test]
    fn s2_res_chain() {
        let mut ctx = s2_res_ctx(10, 10, chain_stack());
        let initial = s2_res_initial(10, &[(1, 1.0), (2, 1.0)], 1.0);
        let mut ws = HighsConflictResolution::default();
        ws.ensure_dim(10);

        let n = analyze_resolution(&mut ctx, &mut ws, &initial, 2);

        assert_eq!(n, 1);
        assert_eq!(ws.conflictrows.len(), 1);
        let row = &ws.conflictrows[0];
        assert_eq!(row.inds, vec![0]);
        assert_eq!(row.vals[0], -1.0);
        assert_eq!(row.lhs, 0.0);
        assert_eq!(ws.nressteps, 2);
        assert_eq!(ws.nfuips, 1);
    }

    #[test]
    fn s2_res_chain_other_order() {
        let stack = vec![
            S2TestChg {
                col: 0,
                is_upper: false,
                newbound: 1.0,
                depth: 1,
                branching: true,
                reason: None,
            },
            S2TestChg {
                col: 1,
                is_upper: true,
                newbound: 0.0,
                depth: 1,
                branching: false,
                reason: Some((vec![(0, -1.0), (1, -1.0)], -1.0)),
            },
            S2TestChg {
                col: 2,
                is_upper: false,
                newbound: 1.0,
                depth: 1,
                branching: false,
                reason: Some((vec![(1, 1.0), (2, 1.0)], 1.0)),
            },
        ];
        let mut ctx = s2_res_ctx(10, 10, stack);
        let initial = s2_res_initial(10, &[(0, -1.0), (2, -1.0)], -1.0);
        let mut ws = HighsConflictResolution::default();
        ws.ensure_dim(10);

        let n = analyze_resolution(&mut ctx, &mut ws, &initial, 2);

        assert_eq!(n, 1);
        assert_eq!(ws.conflictrows.len(), 1);
        let row = &ws.conflictrows[0];
        assert_eq!(row.inds, vec![0]);
        assert_eq!(row.vals[0], -1.0);
        assert_eq!(row.lhs, 0.0);
        assert_eq!(ws.nressteps, 2);
    }

    #[test]
    fn s2_res_all_branchings() {
        let stack = vec![
            S2TestChg {
                col: 0,
                is_upper: false,
                newbound: 3.0,
                depth: 1,
                branching: true,
                reason: None,
            },
            S2TestChg {
                col: 1,
                is_upper: false,
                newbound: 2.0,
                depth: 1,
                branching: true,
                reason: None,
            },
        ];
        let mut ctx = s2_res_ctx(10, 10, stack);
        ctx.gub[0] = 10.0;
        ctx.gub[1] = 10.0;
        let initial = s2_res_initial(10, &[(0, -1.0), (1, -1.0)], -4.0);
        let mut ws = HighsConflictResolution::default();
        ws.ensure_dim(10);

        let n = analyze_resolution(&mut ctx, &mut ws, &initial, 1);

        assert_eq!(n, 0);
        assert!(ws.conflictrows.is_empty());
    }

    #[test]
    fn s2_res_maxsize_small() {
        let mut ctx = s2_res_ctx(4, 4, chain_stack());
        let initial = s2_res_initial(4, &[(1, 1.0), (2, 1.0)], 1.0);
        let mut ws = HighsConflictResolution::default();
        ws.ensure_dim(4);

        let n = analyze_resolution(&mut ctx, &mut ws, &initial, 2);

        assert_eq!(n, 0);
        assert!(ws.conflictrows.is_empty());
        assert!(ws.stats.nreslongconfs >= 1);
    }
}
