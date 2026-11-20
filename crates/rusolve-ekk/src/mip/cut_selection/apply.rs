//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! candidate collection, metric computation, hybrid selection, and commit
//! for the hybrid cut selector. Mirrors SCIPsepastoreApplyCuts + SCIPcutselsSelect.

use crate::mip::cut_pool::HighsCutSet;
use crate::mip::domain::HighsDomain;
use crate::mip::mip_solver::HighsMipSolver;
use crate::prelude::*;

use super::S4_INFINITY;
use super::S4_MAXCUTS;
use super::S4_MAXCUTSROOT;
use super::S4Cand;
use super::S4RandNumGen;
use super::hybrid::s4_cutsel_select_hybrid;
use super::row_metrics::*;

/// SCIP 11.0.0 (32dd425): sepastore.c::SCIPsepastoreApplyCuts (lines 934-1058)
/// HiGHS: HighsSeparation.cpp::separate (lines 154-160)
pub(crate) fn s4_separate_cutpools(
    ms: &mut HighsMipSolver,
    sol: &[f64],
    propdomain: &HighsDomain,
    cutset: &mut HighsCutSet,
    rng: &mut S4RandNumGen,
    feastol: f64,
    cutpool: usize,
    root: bool,
) {
    // 1. Collect candidates from the current pool and (if not pool 0) the global pool.
    let w: Vec<(f64, HighsInt)> = crate::detach!(ms.mipdata_.cutpools[cutpool], |cp| cp
        .s4_collect_candidates(sol, propdomain, feastol, false));

    let g: Vec<(f64, HighsInt)> = if cutpool != 0 {
        crate::detach!(ms.mipdata_.cutpools[0], |cp| cp
            .s4_collect_candidates(sol, propdomain, feastol, true))
    } else {
        Vec::new()
    };

    // 2. Forced cuts already in cutset (sepastore.c:529 forced list; expected 0, Q11).
    let nforced = cutset.cutindices.len();

    // 3. Compute metrics for all candidates (forced + w + g).
    let cost = ms.col_cost();
    let objsqrnorm = s4_obj_sqrnorm(cost);
    let has_incumbent = !ms.mipdata_.incumbent.is_empty();
    let dir: Vec<f64> = if has_incumbent {
        s4_sol_direction(&ms.mipdata_.incumbent, sol)
    } else {
        Vec::new()
    };

    let mut cands: Vec<S4Cand> = Vec::new();
    let mut rows: Vec<(HighsInt, HighsInt)> = Vec::new();

    // Forced cuts: pool and row from cutset.
    for k in 0..nforced {
        let p = cutset.cutpools[k];
        let r = cutset.cutindices[k];
        let cand = compute_cand(ms, p, r, sol, cost, objsqrnorm, &dir, has_incumbent);
        cands.push(cand);
        rows.push((p, r));
    }

    // Candidates from w (pool = cutpool).
    for &(_, row) in &w {
        let cand = compute_cand(
            ms,
            cutpool as HighsInt,
            row,
            sol,
            cost,
            objsqrnorm,
            &dir,
            has_incumbent,
        );
        cands.push(cand);
        rows.push((cutpool as HighsInt, row));
    }

    // Candidates from g (pool = 0).
    for &(_, row) in &g {
        let cand = compute_cand(ms, 0, row, sol, cost, objsqrnorm, &dir, has_incumbent);
        cands.push(cand);
        rows.push((0, row));
    }

    // 4. Selection.
    let maxsepacuts: usize = if root {
        S4_MAXCUTSROOT as usize
    } else {
        S4_MAXCUTS as usize
    };
    let total = cands.len();

    let selected: Vec<usize> = if nforced >= std::cmp::min(total, maxsepacuts) {
        // All forced, no selection needed.
        (0..nforced).collect()
    } else {
        // SCIP cutsel.c:189-190: maxnselectedcuts -= nforcedcuts; MIN(maxnselectedcuts, ncuts).
        // ncuts here is the total count (forced included), as in SCIP. A bound above
        // the non-forced count never binds: the greedy loop ends when no cut is left.
        let maxsel = std::cmp::min(maxsepacuts - nforced, total);

        // Build the parallelism closure. It borrows ms immutably for pool data.
        // rng is separate, so no conflict.
        let par = |a: usize, b: usize| -> f64 {
            let (pa, ra) = rows[a];
            let (pb, rb) = rows[b];
            let pool_a = &ms.mipdata_.cutpools[pa as usize];
            let pool_b = &ms.mipdata_.cutpools[pb as usize];
            let sa = pool_a.matrix_.get_row_start(ra);
            let ea = pool_a.matrix_.get_row_end(ra);
            let sb = pool_b.matrix_.get_row_start(rb);
            let eb = pool_b.matrix_.get_row_end(rb);
            let idx_a = &pool_a.matrix_.get_a_rindex()[sa as usize..ea as usize];
            let val_a = &pool_a.matrix_.get_a_rvalue()[sa as usize..ea as usize];
            let idx_b = &pool_b.matrix_.get_a_rindex()[sb as usize..eb as usize];
            let val_b = &pool_b.matrix_.get_a_rvalue()[sb as usize..eb as usize];
            s4_row_parallelism(idx_a, val_a, idx_b, val_b)
        };

        let sel = s4_cutsel_select_hybrid(
            &cands,
            nforced,
            root,
            Some(rng),
            has_incumbent,
            maxsel,
            &par,
        );

        // Prepend forced indices.
        let mut result: Vec<usize> = (0..nforced).collect();
        result.extend(sel);
        result
    };

    // 5. Apply loop: forced cuts are already in cutset and are not committed again;
    // non-forced cut i is kept iff efficacy > feastol (SCIPsetIsFeasPositive).
    let kept: Vec<usize> = selected
        .into_iter()
        .filter(|&i| i >= nforced && cands[i].efficacy > feastol)
        .collect();

    // 6. Commit: group kept indices into maximal consecutive runs with the same pool.
    if kept.is_empty() {
        return;
    }

    let mut run_start = 0;
    let mut run_pool = rows[kept[0]].0;

    for idx in 1..=kept.len() {
        let end_run = idx == kept.len() || rows[kept[idx]].0 != run_pool;
        if end_run {
            let rows_of_run: Vec<HighsInt> =
                kept[run_start..idx].iter().map(|&i| rows[i].1).collect();
            let ts = run_pool != cutpool as HighsInt;
            crate::detach!(ms.mipdata_.cutpools[run_pool as usize], |cp| {
                cp.s4_commit_selected(cutset, &rows_of_run, ts)
            });
            if idx < kept.len() {
                run_start = idx;
                run_pool = rows[kept[idx]].0;
            }
        }
    }
}

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::scoring (lines 102-188) row metrics, precomputed per D4
/// (SCIProwGetLPEfficacy lp.c:7054-7092, SCIProwGetObjParallelism lp.c:8046-8075,
/// SCIPgetRowNumIntCols/SCIProwGetNNonz, SCIProwGetLPSolCutoffDistance lp.c:6997-7051).
fn compute_cand(
    ms: &HighsMipSolver,
    p: HighsInt,
    r: HighsInt,
    sol: &[f64],
    cost: &[f64],
    objsqrnorm: f64,
    dir: &[f64],
    has_incumbent: bool,
) -> S4Cand {
    let pool = &ms.mipdata_.cutpools[p as usize];
    let start = pool.matrix_.get_row_start(r) as usize;
    let end = pool.matrix_.get_row_end(r) as usize;
    let idx = &pool.matrix_.get_a_rindex()[start..end];
    let val = &pool.matrix_.get_a_rvalue()[start..end];
    let rhs = pool.rhs_[r as usize];

    let act = s4_row_activity(idx, val, sol);
    let feas = s4_row_feasibility(-S4_INFINITY, rhs, act);
    let norm = s4_row_norm(val);
    let efficacy = s4_row_efficacy(feas, norm);

    let nnz = (end - start) as HighsInt;
    debug_assert!(nnz > 1);

    let objpar = s4_row_obj_parallelism(idx, val, cost, objsqrnorm);
    let numint = s4_row_num_int_cols(idx, ms.variable_type());

    let cutoffdist = if has_incumbent {
        s4_row_cutoff_distance(idx, val, rhs, dir, act)
    } else {
        0.0
    };

    let in_global = p == 0;

    S4Cand {
        pool: p,
        row: r,
        efficacy,
        objpar,
        numintcols: numint as f64,
        nnz: nnz as f64,
        cutoffdist,
        in_global,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::Slot;
    use crate::mip::cut_pool::HighsCutPool;

    fn build_fixture(ncols: usize) -> HighsMipSolver {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = ncols as HighsInt;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; ncols];
        lp.col_lower_ = vec![0.0; ncols];
        lp.col_upper_ = vec![1.0; ncols];
        lp.integrality_ = vec![HighsVarType::kInteger; ncols];
        lp.a_matrix_.num_col_ = ncols as HighsInt;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; ncols + 1];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms
    }

    fn add(
        pool: &mut HighsCutPool,
        ms: &mut HighsMipSolver,
        idx: &[HighsInt],
        val: &[f64],
        rhs: f64,
    ) -> HighsInt {
        let mut i = idx.to_vec();
        let mut v = val.to_vec();
        let r = pool.add_cut(
            ms,
            &mut i,
            &mut v,
            idx.len() as HighsInt,
            rhs,
            false,
            false,
            false,
            false,
        );
        assert!(r >= 0);
        r
    }

    /// 150 cuts on disjoint column pairs: parallelism 0, so only maxsepacuts binds.
    fn disjoint_cuts_run(root: bool) -> usize {
        let ncuts = 150;
        let mut ms = build_fixture(2 * ncuts);
        let mut pool0 = HighsCutPool::new(2 * ncuts as HighsInt, 10, 10000, 0);
        for k in 0..ncuts as HighsInt {
            add(
                &mut pool0,
                &mut ms,
                &[2 * k, 2 * k + 1],
                &[-1.0, -1.0],
                -0.8,
            );
        }
        ms.mipdata_.cutpools.push(Slot::new(pool0));
        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 2 * ncuts];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 0, root);
        cutset.cutindices.len()
    }

    #[test]
    fn s4_apply_maxcuts_root_vs_nonroot() {
        assert_eq!(disjoint_cuts_run(true), 150);
        assert_eq!(disjoint_cuts_run(false), 100);
    }

    #[test]
    fn s4_apply_feaspositive_drops_small_efficacy() {
        let mut ms = build_fixture(4);
        let mut pool0 = HighsCutPool::new(4, 10, 100, 0);
        // V: -x0 + x1 <= -0.5. HiGHS score 0.5, SCIP efficacy 0.5/sqrt(2).
        let v = add(&mut pool0, &mut ms, &[0, 1], &[-1.0, 1.0], -0.5);
        // T: -x2 + 1e6 x3 <= -0.48. HiGHS score 0.48 (x3 at its lower bound is not
        // active), so the collect step keeps it, but SCIP efficacy ~4.8e-7 <= feastol.
        let _t = add(&mut pool0, &mut ms, &[2, 3], &[-1.0, 1e6], -0.48);
        ms.mipdata_.cutpools.push(Slot::new(pool0));

        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 4];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 0, true);

        // Both were selected (two RNG draws, both candidates scored) but only V is applied.
        let mut expected_rng = S4RandNumGen::new(0x5EED);
        expected_rng.get_rand();
        expected_rng.get_rand();
        assert_eq!(rng.seed, expected_rng.seed);
        assert_eq!(cutset.cutindices, vec![v]);
        assert_eq!(cutset.cutpools, vec![0]);
    }

    #[test]
    fn s4_apply_cutset_rows_equal_pool_rows() {
        let mut ms = build_fixture(4);
        let mut pool0 = HighsCutPool::new(4, 10, 100, 0);
        add(&mut pool0, &mut ms, &[0, 1], &[-1.0, -2.0], -0.8);
        add(&mut pool0, &mut ms, &[2, 3], &[-3.0, -4.0], -1.8);
        ms.mipdata_.cutpools.push(Slot::new(pool0));
        ms.mipdata_
            .cutpools
            .push(Slot::new(HighsCutPool::new(4, 10, 100, 1)));

        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 4];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 0, true);

        assert_eq!(cutset.cutindices.len(), 2);
        for k in 0..cutset.cutindices.len() {
            let r = cutset.cutindices[k];
            let p = cutset.cutpools[k];
            assert_eq!(p, 0);
            let pool = &ms.mipdata_.cutpools[p as usize];
            let start = pool.matrix_.get_row_start(r) as usize;
            let end = pool.matrix_.get_row_end(r) as usize;
            let cs = cutset.a_rstart_[k] as usize;
            let ce = cutset.a_rstart_[k + 1] as usize;
            assert_eq!(
                &cutset.a_rindex_[cs..ce],
                &pool.matrix_.get_a_rindex()[start..end]
            );
            assert_eq!(
                &cutset.a_rvalue_[cs..ce],
                &pool.matrix_.get_a_rvalue()[start..end]
            );
            assert_eq!(cutset.upper_[k], pool.rhs_[r as usize]);
            assert_eq!(pool.ages_[r as usize], -1);
        }
    }

    #[test]
    fn s4_apply_cutpool_zero_collects_pool0_only() {
        let mut ms = build_fixture(4);
        let mut pool0 = HighsCutPool::new(4, 10, 100, 0);
        let r0 = add(&mut pool0, &mut ms, &[0, 1], &[-1.0, -1.0], -0.8);
        ms.mipdata_.cutpools.push(Slot::new(pool0));

        let mut pool1 = HighsCutPool::new(4, 10, 100, 1);
        let r1 = add(&mut pool1, &mut ms, &[2, 3], &[-1.0, -1.0], -0.8);
        ms.mipdata_.cutpools.push(Slot::new(pool1));
        let age1 = ms.mipdata_.cutpools[1].ages_[r1 as usize];

        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 4];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 0, true);

        assert_eq!(cutset.cutindices, vec![r0]);
        assert_eq!(cutset.cutpools, vec![0]);
        // pool 1 was not touched: its cut is neither aged nor in the LP.
        assert_eq!(ms.mipdata_.cutpools[1].ages_[r1 as usize], age1);
        assert_eq!(ms.mipdata_.cutpools[1].num_lp_cuts, 0);
    }

    #[test]
    fn s4_apply_rng_draws_ncands() {
        let mut ms = build_fixture(4);
        let mut pool0 = HighsCutPool::new(4, 10, 100, 0);
        add(&mut pool0, &mut ms, &[0, 1], &[-1.0, -1.0], -0.8);
        add(&mut pool0, &mut ms, &[1, 2], &[-1.0, -1.0], -0.8);
        add(&mut pool0, &mut ms, &[2, 3], &[-1.0, -1.0], -0.8);
        ms.mipdata_.cutpools.push(Slot::new(pool0));
        ms.mipdata_
            .cutpools
            .push(Slot::new(HighsCutPool::new(4, 10, 100, 1)));

        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 4];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        let ncands = 3;
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 0, true);

        let mut expected_rng = S4RandNumGen::new(0x5EED);
        for _ in 0..ncands {
            expected_rng.get_rand();
        }
        assert_eq!(rng.seed, expected_rng.seed);
        assert_eq!(rng.xor_seed, expected_rng.xor_seed);
        assert_eq!(rng.mwc_seed, expected_rng.mwc_seed);
        assert_eq!(rng.cst_seed, expected_rng.cst_seed);

        // No candidates: no selection, no draw.
        let mut cutset2 = HighsCutSet::default();
        let before = rng.clone();
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset2, &mut rng, 1e-6, 0, true);
        assert!(cutset2.cutindices.is_empty());
        assert_eq!(rng.seed, before.seed);
        assert_eq!(rng.mwc_seed, before.mwc_seed);
    }

    #[test]
    fn s4_apply_commit_order_interleaved() {
        let mut ms = build_fixture(8);
        // Disjoint supports, each cut -x_j + w x_k <= -1. At sol 0, x_k sits at its lower
        // bound, so every HiGHS score is 1 (no pool truncation), while the SCIP efficacy is
        // 1/sqrt(1+w^2): A > B > C > D with A, C in pool 1 and B, D in pool 0, so the
        // selected order alternates pools: four single-cut runs.
        let mut pool0 = HighsCutPool::new(8, 10, 100, 0);
        let rb = add(&mut pool0, &mut ms, &[2, 3], &[-1.0, 1.0], -1.0);
        let rd = add(&mut pool0, &mut ms, &[6, 7], &[-1.0, 3.0], -1.0);
        ms.mipdata_.cutpools.push(Slot::new(pool0));

        let mut pool1 = HighsCutPool::new(8, 10, 100, 1);
        let ra = add(&mut pool1, &mut ms, &[0, 1], &[-1.0, 0.5], -1.0);
        let rc = add(&mut pool1, &mut ms, &[4, 5], &[-1.0, 2.0], -1.0);
        ms.mipdata_.cutpools.push(Slot::new(pool1));

        let age_b = ms.mipdata_.cutpools[0].ages_[rb as usize];
        let dom = HighsDomain::new(&ms);
        let sol = vec![0.0; 8];
        let mut cutset = HighsCutSet::default();
        let mut rng = S4RandNumGen::new(0x5EED);
        s4_separate_cutpools(&mut ms, &sol, &dom, &mut cutset, &mut rng, 1e-6, 1, true);

        assert_eq!(cutset.cutpools, vec![1, 0, 1, 0]);
        assert_eq!(cutset.cutindices, vec![ra, rb, rc, rd]);
        // pool 1 (thread_safe = false) books its cuts as LP cuts; pool 0 (thread_safe = true) does not.
        assert_eq!(ms.mipdata_.cutpools[1].num_lp_cuts, 2);
        assert_eq!(ms.mipdata_.cutpools[1].ages_[ra as usize], -1);
        assert_eq!(ms.mipdata_.cutpools[1].ages_[rc as usize], -1);
        assert_eq!(ms.mipdata_.cutpools[0].num_lp_cuts, 0);
        assert_eq!(ms.mipdata_.cutpools[0].ages_[rb as usize], age_b);
    }
}
