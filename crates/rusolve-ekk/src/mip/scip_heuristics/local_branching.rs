//! Translated from SCIP (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! `mip_local_branching`: local branching on the port stack.
//!
//! * M. Fischetti, A. Lodi, "Local branching", Math. Programming 98 (2003).
//! * SCIP `src/scip/heur_localbranching.c` (63417ae): `heurExecLocalbranching` (gates, the node
//!   limit `nodesquot * nnodes`, rewarded by the success rate, `-100` nodes per call, `+ nodesofs`),
//!   `addLocalbranchingConstraintAndObjcutoff` (Hamming ball of radius `curneighborhoodsize`
//!   around the incumbent on the binaries, objective cutoff `(1 - minimprove) ub + minimprove lb`),
//!   `setupAndSolveSubscipLocalbranching` (radius and node-minimum adaptation on the sub-MIP
//!   status). Defaults: radius 18, nodesofs 1000, maxnodes 10000, minimprove 0.01, minnodes 1000,
//!   nodesquot 0.05, waitingnodes 200, bestsollimit 3.
//!
//! Divergences: the sub-MIP is the sub-MIP copy (`s6_solve_sub_mip_ext`); the LP-iteration
//! event handler and `copycuts` have no counterpart; SCIP's "stall nodes" are
//! `max(10, nodes / 10)` as in `SCIPsetCommonSubscipParams`. The status cases map as: found or
//! solved to optimality -> wait for a new solution; infeasible -> grow the radius; out of nodes
//! without a solution -> shrink it and double the node minimum. No randomness.

use super::proximity::{
    append_row, basis_with_extra_row, binary_cols, incumbent_of, rounded_lower_bound,
};
use super::sub_mip_s6::S6SubMipExt;
use super::*;
use crate::mip::mip_solver::HighsMipSolver;
use crate::mip::mip_solver_data::kSolutionSourceLocalBranching;
use crate::work::WorkSnapshot;

/// SCIP heur_localbranching.c DEFAULT_NEIGHBORHOODSIZE
pub(crate) const LB_NEIGHBORHOODSIZE: i64 = 18;
/// SCIP heur_localbranching.c DEFAULT_NODESOFS
pub(crate) const LB_NODESOFS: i64 = 1000;
/// SCIP heur_localbranching.c DEFAULT_MAXNODES
pub(crate) const LB_MAXNODES: i64 = 10000;
/// SCIP heur_localbranching.c DEFAULT_MINIMPROVE
pub(crate) const LB_MINIMPROVE: f64 = 0.01;
/// SCIP heur_localbranching.c DEFAULT_MINNODES
pub(crate) const LB_MINNODES: i64 = 1000;
/// SCIP heur_localbranching.c DEFAULT_NODESQUOT
pub(crate) const LB_NODESQUOT: f64 = 0.05;
/// SCIP heur_localbranching.c DEFAULT_NWAITINGNODES
pub(crate) const LB_WAITINGNODES: i64 = 200;
/// SCIP heur_localbranching.c DEFAULT_BESTSOLLIMIT
pub(crate) const LB_BESTSOLLIMIT: HighsInt = 3;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum LbCallStatus {
    Execute,
    WaitForNewSol,
}

#[derive(Clone, Debug)]
pub(crate) struct LocalBranchingState {
    pub curneighborhoodsize: i64,
    pub emptyneighborhoodsize: i64,
    pub curminnodes: i64,
    pub usednodes: i64,
    pub callstatus: LbCallStatus,
    /// solution index of the incumbent the radius state belongs to
    pub lastsol: Option<i64>,
    /// `SCIPheurGetNCalls` / `SCIPheurGetNBestSolsFound`
    pub ncalls: i64,
    pub nbestsolsfound: i64,
}

impl Default for LocalBranchingState {
    fn default() -> Self {
        Self {
            curneighborhoodsize: LB_NEIGHBORHOODSIZE,
            emptyneighborhoodsize: 0,
            curminnodes: LB_MINNODES,
            usednodes: 0,
            callstatus: LbCallStatus::Execute,
            lastsol: None,
            ncalls: 0,
            nbestsolsfound: 0,
        }
    }
}

/// SCIP `heurExecLocalbranching`.
pub(crate) fn local_branching_step(ms: &mut HighsMipSolver) {
    if ms.submip || !ms.options_mip_.mip_local_branching || ms.terminate() {
        return;
    }
    if ms.mipdata_.nodequeue.empty() || ms.mipdata_.local_branching.is_none() {
        return;
    }
    let bins = binary_cols(ms);
    // enough binary variables that a local branching constraint makes sense
    if (bins.len() as i64) < 2 * LB_NEIGHBORHOODSIZE {
        return;
    }
    let Some(inc) = incumbent_of(ms) else {
        return;
    };
    let nnodes = ms.mipdata_.num_nodes;
    // enough nodes since the last incumbent
    if nnodes - inc.nodenum < LB_WAITINGNODES {
        return;
    }
    let mut st = ms.mipdata_.local_branching.take().expect("checked above");
    local_branching_inner(ms, &mut st, &bins, &inc, nnodes);
    ms.mipdata_.local_branching = Some(st);
}

fn local_branching_inner(
    ms: &mut HighsMipSolver,
    st: &mut LocalBranchingState,
    bins: &[usize],
    inc: &super::proximity::Incumbent,
    nnodes: i64,
) {
    // reset neighbourhood and minnodes if a new solution was found
    if st.lastsol != Some(inc.index) {
        st.curneighborhoodsize = LB_NEIGHBORHOODSIZE;
        st.curminnodes = LB_MINNODES;
        st.emptyneighborhoodsize = 0;
        st.callstatus = LbCallStatus::Execute;
        st.lastsol = Some(inc.index);
    }
    if st.callstatus == LbCallStatus::WaitForNewSol {
        return;
    }
    // maximal number of branching nodes until the heuristic is aborted
    let mut maxnnodes = (LB_NODESQUOT * nnodes as f64) as i64;
    // reward local branching if it succeeded often
    maxnnodes = (maxnnodes as f64
        * (1.0 + 2.0 * (st.nbestsolsfound as f64 + 1.0) / (st.ncalls as f64 + 1.0)))
        as i64;
    // the setup costs of the sub-MIP count as 100 nodes
    maxnnodes -= 100 * st.ncalls;
    maxnnodes += LB_NODESOFS;
    let nsubnodes = (maxnnodes - st.usednodes).min(LB_MAXNODES);
    if nsubnodes < st.curminnodes {
        return;
    }

    // the local branching row over the binaries and the objective cutoff
    let feastol = ms.mipdata_.feastol;
    let mut lhs = (st.emptyneighborhoodsize + 1) as f64;
    let mut rhs = st.curneighborhoodsize as f64;
    let mut idx: Vec<HighsInt> = Vec::with_capacity(bins.len());
    let mut val: Vec<f64> = Vec::with_capacity(bins.len());
    for &j in bins {
        if inc.x[j] > 1.0 - feastol {
            val.push(-1.0);
            rhs -= 1.0;
            lhs -= 1.0;
        } else {
            val.push(1.0);
        }
        idx.push(j as HighsInt);
    }
    let ub = ms.mipdata_.upper_bound;
    let lb = rounded_lower_bound(ms, ms.mipdata_.lower_bound);
    let upperbound = ub - S6_SUMEPSILON;
    let mut cutoff = if !s6_is_inf(lb.abs()) {
        (1.0 - LB_MINIMPROVE) * ub + LB_MINIMPROVE * lb
    } else if ub >= 0.0 {
        (1.0 - LB_MINIMPROVE) * ub
    } else {
        (1.0 + LB_MINIMPROVE) * ub
    };
    cutoff = cpp_min(upperbound, cutoff);

    let mut lp: HighsLp = mip_model!(ms).clone();
    let num_row = lp.num_row_ as usize;
    append_row(&mut lp, "localbranching", idx, val, lhs, rhs);
    let basis = basis_with_extra_row(ms, num_row);
    let d = &ms.mipdata_.domains[0];
    let (cl, cu) = (d.col_lower_.clone(), d.col_upper_.clone());

    let old_upper = ms.mipdata_.upper_bound;
    let w0 = WorkSnapshot::now();
    let (_ok, subnodes, status) = detach!(ms.mipdata_.heuristics, |h| h.s6_solve_sub_mip_ext(
        ms,
        0,
        lp,
        basis,
        0.0,
        cl,
        cu,
        kHighsIInf,
        nsubnodes as HighsInt,
        cpp_max(10, nsubnodes / 10) as HighsInt,
        cutoff,
        kSolutionSourceLocalBranching,
        S6SubMipExt {
            own_objective: false,
            max_improving_sols: Some(LB_BESTSOLLIMIT),
            ..Default::default()
        },
    ));
    let work = WorkSnapshot::now().since(&w0).work();
    let improved = ms.mipdata_.upper_bound < old_upper;
    st.usednodes += subnodes;
    st.ncalls += 1;
    if improved {
        st.nbestsolsfound += 1;
    }
    let cur = st.curneighborhoodsize;
    match status {
        HighsModelStatus::kInfeasible => {
            // the neighbourhood is empty: grow it
            st.emptyneighborhoodsize = cur;
            st.curneighborhoodsize = cur + cur / 2;
            st.curneighborhoodsize = cpp_max(st.curneighborhoodsize, st.emptyneighborhoodsize + 2);
            st.callstatus = LbCallStatus::Execute;
        }
        HighsModelStatus::kSolutionLimit | HighsModelStatus::kIterationLimit if !improved => {
            // node limit without a solution: shrink the neighbourhood, double the node minimum
            st.callstatus = LbCallStatus::Execute;
            st.curneighborhoodsize = (st.emptyneighborhoodsize + cur) / 2;
            st.curminnodes *= 2;
            if st.curneighborhoodsize <= st.emptyneighborhoodsize {
                st.callstatus = LbCallStatus::WaitForNewSol;
            }
        }
        // solved, improved, or anything unexpected: wait for a new solution
        _ => st.callstatus = LbCallStatus::WaitForNewSol,
    }
    let s = &mut ms.stats_;
    s.lb_calls += 1;
    s.lb_nodes += subnodes;
    s.lb_work += work;
    if improved {
        s.lb_found += 1;
    }
    crate::io::highs_io::highs_log_dev(
        &ms.options_mip_.log_options,
        HighsLogType::kInfo,
        &format!(
            "LocalBranching radius {} nodes {} work {} improved {}\n",
            cur, subnodes, work, improved
        ),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn state_starts_at_the_scip_defaults() {
        let st = LocalBranchingState::default();
        assert_eq!(st.curneighborhoodsize, 18);
        assert_eq!(st.curminnodes, 1000);
        assert_eq!(st.callstatus, LbCallStatus::Execute);
    }
}
