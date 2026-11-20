//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP 11.0.0 heur_crossover.c: crossover primal heuristic.

use super::*;
use crate::mip::mip_solver::HighsMipSolver;

/// SCIP 11.0.0 (32dd425): heur_crossover.c::calculateHashKey (lines 201-218)
pub(crate) fn s6_calculate_hash_key(sel: &[i32]) -> u32 {
    let mut hashkey: u32 = 1;
    for i in 0..sel.len() {
        hashkey = hashkey.wrapping_mul((sel[i] as u32).wrapping_add(1));
    }
    for i in 0..sel.len() {
        hashkey = hashkey.wrapping_add(sel[i] as u32);
    }
    hashkey
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::sortArray (lines 221-243)
pub(crate) fn s6_sort_array(a: &mut [i32]) {
    let size = a.len();
    for i in 1..size {
        let tmp = a[i];
        let mut j: i32 = i as i32 - 1;
        while j >= 0 && a[j as usize] > tmp {
            a[(j + 1) as usize] = a[j as usize];
            j -= 1;
        }
        a[(j + 1) as usize] = tmp;
    }
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::createSolTuple (lines 247-270)
pub(crate) fn s6_create_sol_tuple(
    selection: &mut [i32],
    lasttuple: &mut Option<Vec<i32>>,
) -> Vec<i32> {
    let tuple: Vec<i32> = selection.to_vec();
    // SCIP quirk (heur_crossover.c:247-270, :774-788)
    s6_sort_array(selection);
    *lasttuple = Some(tuple.clone());
    tuple
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::solHasNewSource (lines 274-290)
pub(crate) fn s6_sol_has_new_source(
    pool: &S6SolPool,
    selection: &[i32],
    selectionsize: usize,
    newsol: usize,
) -> bool {
    let new_origin = pool.sols[newsol].origin;
    let new_nodenum = pool.sols[newsol].nodenum;
    for i in 0..selectionsize {
        let sel_idx = selection[i] as usize;
        if pool.sols[sel_idx].origin == new_origin && pool.sols[sel_idx].nodenum == new_nodenum {
            return false;
        }
    }
    true
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::selectSolsRandomized (lines 295-357)
pub(crate) fn s6_select_sols_randomized(
    st: &mut S6State,
    selection: &mut [i32],
    nusedsols: usize,
) -> bool {
    debug_assert!(st.xover_rng.is_some());
    let mut lastsol: i32 = st.pool.sols.len() as i32;
    let mut success = false;
    let mut i: i32 = 0;

    while !success && i < S6_XOVER_SELECT_TRIES {
        let mut validtuple = true;
        let mut j: usize = 0;
        while j < nusedsols && validtuple {
            let rng = st.xover_rng.as_mut().unwrap();
            let mut k: i32 = s6_get_int(rng, (nusedsols - j - 1) as i32, lastsol - 1);

            while k >= (nusedsols - j - 1) as i32
                && !s6_sol_has_new_source(&st.pool, selection, j, k as usize)
            {
                k -= 1;
            }

            validtuple = k >= (nusedsols - j - 1) as i32;
            selection[j] = k;
            lastsol = k;
            j += 1;
        }

        if validtuple {
            let elem = s6_create_sol_tuple(&mut selection[..nusedsols], &mut st.xover_lasttuple);
            if !st.xover_hashtable.contains(&elem) {
                st.xover_hashtable.insert(elem);
                success = true;
            }
        }
        i += 1;
    }

    success
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::fixVariables (lines 362-436)
pub(crate) fn s6_fix_variables(
    st: &S6State,
    integrality: &[HighsVarType],
    col_lower: &[f64],
    col_upper: &[f64],
    selection: &[i32],
    nusedsols: usize,
    fixed: &mut Vec<(usize, f64)>,
) -> bool {
    fixed.clear();

    let nint: usize = integrality
        .iter()
        .filter(|&&t| t == HighsVarType::kInteger)
        .count();

    for col in 0..integrality.len() {
        if integrality[col] != HighsVarType::kInteger {
            continue;
        }

        let mut fixable = true;
        let solval: f64 = st.pool.sols[selection[0] as usize].x[col];

        for j in 1..nusedsols {
            let varsolval: f64 = st.pool.sols[selection[j] as usize].x[col];
            if (solval - varsolval).abs() > 0.5 {
                fixable = false;
                break;
            }
        }

        fixable = fixable && col_lower[col] <= solval && solval <= col_upper[col];

        if fixable {
            fixed.push((col, solval));
        }
    }

    let fixingrate: f64 = fixed.len() as f64 / (cpp_max(nint as i64, 1i64) as usize) as f64;

    fixed.len() < nint && fixingrate >= S6_XOVER_MINFIXINGRATE
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::determineVariableFixings (lines 440-518)
pub(crate) fn s6_determine_variable_fixings(
    st: &mut S6State,
    integrality: &[HighsVarType],
    col_lower: &[f64],
    col_upper: &[f64],
    selection: &mut [i32],
    nusedsols: usize,
    fixed: &mut Vec<(usize, f64)>,
) -> bool {
    let nsols: usize = st.pool.sols.len();
    let randomization: bool = S6_XOVER_RANDOMIZATION;

    debug_assert!(nusedsols > 1);
    debug_assert!(nsols >= nusedsols);

    let mut success: bool;

    if !randomization
        || nsols == nusedsols
        || st.xover_prevlastsol != Some(st.pool.sols[nusedsols - 1].index)
    {
        for i in 0..nusedsols {
            selection[i] = i as i32;
        }
        let elem = s6_create_sol_tuple(&mut selection[..nusedsols], &mut st.xover_lasttuple);

        let solheur_origin = st.pool.sols[0].origin;
        let solnodenum = st.pool.sols[0].nodenum;
        let mut allsame = true;

        for i in 1..nusedsols {
            if st.pool.sols[i].origin != solheur_origin || st.pool.sols[i].nodenum != solnodenum {
                allsame = false;
            }
        }
        success = !allsame && !st.xover_hashtable.contains(&elem);

        if !st.xover_hashtable.contains(&elem) {
            st.xover_hashtable.insert(elem);
        }

        if !success && randomization && nsols > nusedsols {
            success = s6_select_sols_randomized(st, selection, nusedsols);
        }
    } else {
        success = s6_select_sols_randomized(st, selection, nusedsols);
    }

    if !success {
        return false;
    }

    s6_fix_variables(
        st,
        integrality,
        col_lower,
        col_upper,
        selection,
        nusedsols,
        fixed,
    )
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::updateFailureStatistic (lines 522-532)
pub(crate) fn s6_update_failure_statistic(st: &mut S6State, nnodes: i64) {
    st.xover_nfailures += 1;
    st.xover_nextnodenumber = if st.xover_nfailures <= S6_XOVER_MAXFAILURES_SHIFT {
        nnodes + 100 * (2i64 << st.xover_nfailures)
    } else {
        i64::MAX
    };
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::setupAndSolveSubscipCrossover (lines 700-720, offset-free)
pub(crate) fn s6_xover_cutoff(upper_bound: f64, lower_bound: f64) -> f64 {
    let upperbound: f64 = upper_bound - S6_SUMEPSILON;
    let cutoff: f64;
    if !s6_is_inf(-lower_bound) {
        cutoff = (1.0 - S6_XOVER_MINIMPROVE) * upper_bound + S6_XOVER_MINIMPROVE * lower_bound;
    } else {
        if upper_bound >= 0.0 {
            cutoff = (1.0 - S6_XOVER_MINIMPROVE) * upper_bound;
        } else {
            cutoff = (1.0 + S6_XOVER_MINIMPROVE) * upper_bound;
        }
    }
    cpp_min(upperbound, cutoff)
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::heurExecCrossover (lines 950-968, arithmetic only)
pub(crate) fn s6_xover_nstallnodes(
    nnodes: i64,
    nbestsolsfound: i64,
    ncalls: i64,
    nodesofs: i64,
    usednodes: i64,
    maxnodes: i64,
) -> i64 {
    let mut nstall: i64 = (S6_XOVER_NODESQUOT * nnodes as f64) as i64;
    nstall = (nstall as f64 * (1.0 + 2.0 * (nbestsolsfound as f64 + 1.0) / (ncalls as f64 + 1.0)))
        as i64;
    nstall -= 100 * ncalls;
    nstall += nodesofs;
    nstall -= usednodes;
    nstall = cpp_min(nstall, maxnodes);
    nstall
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::heurExecCrossover (lines 920-921)
pub(crate) fn s6_xover_enough_sols(nsols: usize, nusedsols: usize) -> bool {
    nsols >= nusedsols
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::setupAndSolveSubscipCrossover (lines 586-811, as mapped in D8)
#[allow(clippy::too_many_arguments)]
pub(crate) fn s6_setup_and_solve_subscip_crossover(
    ms: &mut HighsMipSolver,
    selection: &mut [i32],
    nusedsols: usize,
    fixed: &[(usize, f64)],
    fixingrate: f64,
    nstallnodes: i64,
) -> S6Result {
    // SCIPcopyLargeNeighborhoodSearch, the sub-SCIP parameter settings (catchctrlc, verblevel,
    // limits/bestsol -1, separating off, fast presolving, estimate node selection, inference
    // branching, conflict settings, lp/checkdualfeas) and the LP event handler are not ported:
    // the sub-MIP is HiGHS solveSubMip on *model_ with its own settings.
    // permute is FALSE (heur_crossover.c:81), so SCIPpermuteProb never runs.
    let ub: f64 = ms.mipdata_.upper_bound;
    let lb: f64 = ms.mipdata_.lower_bound;
    let mut lower: Vec<f64> = ms.mipdata_.domains[0].col_lower_.clone();
    let mut upper: Vec<f64> = ms.mipdata_.domains[0].col_upper_.clone();
    let lp: HighsLp = mip_model!(ms).clone();
    let basis: HighsBasis = ms.mipdata_.firstrootbasis.clone();

    debug_assert!(!s6_is_inf(ub));

    let cutoff: f64 = s6_xover_cutoff(ub, lb);

    for &(col, val) in fixed {
        lower[col] = val;
        upper[col] = val;
    }

    let idx_before: i64 = ms.mipdata_.s6.as_mut().unwrap().next_sol_index;
    ms.mipdata_.s6.as_mut().unwrap().cur_origin = S6Origin::Crossover;

    let (_found, subnodes): (bool, i64) = detach!(ms.mipdata_.heuristics, |h| h.s6_solve_sub_mip(
        ms,
        0,
        lp,
        basis,
        fixingrate,
        lower,
        upper,
        kHighsIInf,
        nstallnodes as HighsInt,
        kHighsIInf,
        cutoff
    ));

    let st = ms.mipdata_.s6.as_mut().unwrap();
    st.xover_usednodes += subnodes;
    // SCIPmergeVariableStatistics is not ported: HiGHS solveSubMip has no history transfer.
    st.cur_origin = S6Origin::default();

    let solindex: Option<i64> = st
        .pool
        .sols
        .iter()
        .filter(|s| s.index >= idx_before && s.origin == S6Origin::Crossover)
        .map(|s| s.index)
        .min();

    let mut result: S6Result = S6Result::DidNotFind;

    if let Some(solindex) = solindex {
        result = S6Result::FoundSol;

        for i in 0..nusedsols {
            let tmp: i32 = selection[i];
            selection[i] = solindex as i32;
            let elem = s6_create_sol_tuple(&mut selection[..nusedsols], &mut st.xover_lasttuple);
            st.xover_hashtable.insert(elem);
            selection[i] = tmp;
        }

        if !S6_XOVER_RANDOMIZATION {
            st.xover_prevbestsol = Some(st.pool.sols[0].index);
            st.xover_prevlastsol = Some(st.pool.sols[nusedsols - 1].index);
        }
    }

    let best_index: Option<i64> = st.pool.best().map(|s| s.index);
    if solindex.is_none() || solindex != best_index {
        s6_update_failure_statistic(st, ms.mipdata_.num_nodes);
    }

    result
}

/// SCIP 11.0.0 (32dd425): heur_crossover.c::heurExecCrossover (lines 903-1044)
pub(crate) fn s6_crossover_exec(ms: &mut HighsMipSolver) -> S6Result {
    let nusedsols: usize = S6_XOVER_NUSEDSOLS as usize;
    let nnodes: i64 = ms.mipdata_.num_nodes;

    if ms.mipdata_.s6.is_none() {
        return S6Result::DidNotRun;
    }

    let st = ms.mipdata_.s6.as_mut().unwrap();

    if !s6_xover_enough_sols(st.pool.sols.len(), nusedsols) {
        return S6Result::Delayed;
    }

    if Some(st.pool.sols[nusedsols - 1].index) != st.xover_prevlastsol {
        st.xover_nextnodenumber = nnodes;
        if Some(st.pool.sols[0].index) != st.xover_prevbestsol {
            st.xover_nfailures = 0;
        }
    } else if !S6_XOVER_RANDOMIZATION {
        return S6Result::Delayed;
    }

    if nnodes < st.xover_nextnodenumber {
        return S6Result::Delayed;
    }

    let bestnodenum: i64 = st.pool.sols[0].nodenum;
    // dontwaitatroot is FALSE, so the depth test is always true and is not ported
    if nnodes - bestnodenum < S6_XOVER_NWAITINGNODES {
        return S6Result::Delayed;
    }

    let nstallnodes: i64 = s6_xover_nstallnodes(
        nnodes,
        st.xover_nbestsolsfound,
        st.xover_ncalls,
        S6_XOVER_NODESOFS,
        st.xover_usednodes,
        S6_XOVER_MAXNODES,
    );
    if nstallnodes < S6_XOVER_MINNODES {
        return S6Result::DidNotRun;
    }

    // SCIPcheckCopyLimits and SCIPisStopped are not ported: no decision reads the wall clock
    // (manifest Basics, Determinism), and the main loop checks termination itself.

    let integrality: Vec<HighsVarType> = mip_model!(ms).integrality_.clone();
    let col_lower: Vec<f64> = ms.mipdata_.domains[0].col_lower_.clone();
    let col_upper: Vec<f64> = ms.mipdata_.domains[0].col_upper_.clone();

    let nint: usize = integrality
        .iter()
        .filter(|&&t| t == HighsVarType::kInteger)
        .count();
    if nint == 0 {
        return S6Result::DidNotRun;
    }

    let mut selection: Vec<i32> = vec![0i32; nusedsols];
    let mut fixed: Vec<(usize, f64)> = Vec::new();

    let success: bool = s6_determine_variable_fixings(
        st,
        &integrality,
        &col_lower,
        &col_upper,
        &mut selection,
        nusedsols,
        &mut fixed,
    );

    st.xover_prevbestsol = Some(st.pool.sols[0].index);
    st.xover_prevlastsol = Some(st.pool.sols[nusedsols - 1].index);

    if !success {
        s6_update_failure_statistic(st, nnodes);
        return S6Result::DidNotRun;
    }

    let fixingrate: f64 = fixed.len() as f64 / (cpp_max(nint as i64, 1i64) as usize) as f64;

    let result: S6Result = s6_setup_and_solve_subscip_crossover(
        ms,
        &mut selection,
        nusedsols,
        &fixed,
        fixingrate,
        nstallnodes,
    );

    if result == S6Result::FoundSol || result == S6Result::DidNotFind {
        ms.mipdata_.s6.as_mut().unwrap().xover_ncalls += 1;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s6_xover_create_sol_tuple_unsorted_copy() {
        let mut selection: Vec<i32> = vec![5, 2, 9];
        let mut lasttuple: Option<Vec<i32>> = None;
        let tuple = s6_create_sol_tuple(&mut selection, &mut lasttuple);
        assert_eq!(tuple, vec![5, 2, 9]);
        assert_eq!(selection, vec![2, 5, 9]);
        assert_eq!(lasttuple, Some(vec![5, 2, 9]));
    }

    #[test]
    fn s6_xover_tuple_set_order_sensitive() {
        let mut set: std::collections::HashSet<Vec<i32>> = std::collections::HashSet::new();
        set.insert(vec![1, 2, 3]);
        assert!(set.contains(&vec![1, 2, 3]));
        assert!(!set.contains(&vec![1, 3, 2]));
    }

    #[test]
    fn s6_xover_hash_key_and_sort() {
        assert_eq!(s6_calculate_hash_key(&[0, 1]), 3);
        let mut a: Vec<i32> = vec![3, 1, 2];
        s6_sort_array(&mut a);
        assert_eq!(a, vec![1, 2, 3]);
    }

    #[test]
    fn s6_xover_sol_has_new_source() {
        let mut pool = S6SolPool::default();
        pool.sols.push(S6Sol {
            index: 0,
            obj: 1.0,
            x: vec![0.0],
            origin: S6Origin::Highs(0),
            nodenum: 10,
        });
        pool.sols.push(S6Sol {
            index: 1,
            obj: 2.0,
            x: vec![1.0],
            origin: S6Origin::Highs(1),
            nodenum: 20,
        });
        pool.sols.push(S6Sol {
            index: 2,
            obj: 3.0,
            x: vec![2.0],
            origin: S6Origin::Highs(0),
            nodenum: 10,
        });

        let selection: Vec<i32> = vec![0];
        assert!(!s6_sol_has_new_source(&pool, &selection, 1, 2));

        pool.sols[2].nodenum = 30;
        assert!(s6_sol_has_new_source(&pool, &selection, 1, 2));
    }

    #[test]
    fn s6_xover_select_sols_randomized_seed7() {
        let mut st = S6State::from_flags(true, 0);
        for i in 0..5 {
            st.pool.sols.push(S6Sol {
                index: i as i64,
                obj: i as f64,
                x: vec![0.0],
                origin: S6Origin::Highs(i),
                nodenum: i as i64,
            });
        }
        let nusedsols: usize = 3;
        let mut selection: Vec<i32> = vec![0; nusedsols];
        let success = s6_select_sols_randomized(&mut st, &mut selection, nusedsols);
        assert!(success);

        let mut rng = S4RandNumGen::new(s4_initialize_random_seed_from(S6_XOVER_RANDSEED, 0));
        let k0: i32 = s6_get_int(&mut rng, 2, 4);
        let k1: i32 = s6_get_int(&mut rng, 1, k0 - 1);
        let k2: i32 = s6_get_int(&mut rng, 0, k1 - 1);
        // createSolTuple sorts the caller's selection in place (heur_crossover.c:247-270),
        // so the drawn positions come back in ascending order.
        let mut expected = vec![k0, k1, k2];
        expected.sort();
        assert_eq!(selection, expected);
    }

    #[test]
    fn s6_xover_fix_variables_threshold() {
        let mut st = S6State::from_flags(true, 0);
        let integrality: Vec<HighsVarType> = vec![HighsVarType::kInteger; 4];
        let col_lower: Vec<f64> = vec![0.0; 4];
        let col_upper: Vec<f64> = vec![10.0; 4];
        let nusedsols: usize = 3;

        for i in 0..3 {
            st.pool.sols.push(S6Sol {
                index: i as i64,
                obj: i as f64,
                x: vec![1.0, 2.0, 3.0, 4.0],
                origin: S6Origin::Highs(i),
                nodenum: i as i64,
            });
        }

        let selection: Vec<i32> = vec![0, 1, 2];
        let mut fixed: Vec<(usize, f64)> = Vec::new();

        // (a) 3 identical out of 4 -> true (0.75 >= 0.666)
        st.pool.sols[1].x[3] = 9.0;
        let result = s6_fix_variables(
            &st,
            &integrality,
            &col_lower,
            &col_upper,
            &selection,
            nusedsols,
            &mut fixed,
        );
        assert_eq!(fixed.len(), 3);
        assert!(result);

        // (b) all 4 identical -> false (nfixed == nint)
        st.pool.sols[1].x[3] = 4.0;
        let result = s6_fix_variables(
            &st,
            &integrality,
            &col_lower,
            &col_upper,
            &selection,
            nusedsols,
            &mut fixed,
        );
        assert_eq!(fixed.len(), 4);
        assert!(!result);

        // (c) 2 identical out of 4 -> false (0.5 < 0.666)
        st.pool.sols[1].x[2] = 9.0;
        st.pool.sols[1].x[3] = 9.0;
        let result = s6_fix_variables(
            &st,
            &integrality,
            &col_lower,
            &col_upper,
            &selection,
            nusedsols,
            &mut fixed,
        );
        assert_eq!(fixed.len(), 2);
        assert!(!result);

        // (d) identical column outside bounds is not fixed
        let col_upper_tight: Vec<f64> = vec![10.0, 10.0, 10.0, 2.0];
        st.pool.sols[1].x[2] = 3.0;
        st.pool.sols[1].x[3] = 4.0;
        let result = s6_fix_variables(
            &st,
            &integrality,
            &col_lower,
            &col_upper_tight,
            &selection,
            nusedsols,
            &mut fixed,
        );
        assert_eq!(fixed.len(), 3);
        assert!(!fixed.iter().any(|&(col, _)| col == 3));
    }

    #[test]
    fn s6_xover_update_failure_statistic() {
        let mut st = S6State::from_flags(true, 0);
        s6_update_failure_statistic(&mut st, 1000);
        assert_eq!(st.xover_nfailures, 1);
        assert_eq!(st.xover_nextnodenumber, 1400);

        st.xover_nfailures = 24;
        s6_update_failure_statistic(&mut st, 1000);
        assert_eq!(st.xover_nfailures, 25);
        assert_eq!(st.xover_nextnodenumber, 1000 + 100 * (2i64 << 25));

        s6_update_failure_statistic(&mut st, 1000);
        assert_eq!(st.xover_nfailures, 26);
        assert_eq!(st.xover_nextnodenumber, i64::MAX);
    }

    #[test]
    fn s6_xover_nstallnodes_arithmetic() {
        assert_eq!(s6_xover_nstallnodes(1000, 0, 0, 500, 0, 5000), 800);
        assert!(s6_xover_nstallnodes(100, 0, 0, 500, 1000, 5000) < S6_XOVER_MINNODES);
    }

    #[test]
    fn s6_xover_cutoff_finite_and_infinite_lb() {
        let cut = s6_xover_cutoff(100.0, 0.0);
        assert!((cut - 99.0).abs() < 1e-9);

        let cut = s6_xover_cutoff(100.0, -S6_INFINITY);
        assert!((cut - 99.0).abs() < 1e-9);

        let cut = s6_xover_cutoff(-10.0, -S6_INFINITY);
        assert!((cut - (-10.1)).abs() < 1e-9);
    }

    #[test]
    fn s6_xover_exec_delayed_few_sols() {
        assert!(!s6_xover_enough_sols(2, 3));
        assert!(s6_xover_enough_sols(3, 3));
    }
}
