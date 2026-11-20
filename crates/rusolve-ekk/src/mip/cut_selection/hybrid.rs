//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::row_metrics::s4_max;
use super::*;

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::scoring (lines 82-194)
pub(crate) fn s4_scoring(
    cands: &[S4Cand],
    rng: Option<&mut S4RandNumGen>,
    dircutoffdistweight: f64,
    efficacyweight: f64,
    objparalweight: f64,
    intsupportweight: f64,
    has_incumbent: bool,
    mut scores: Option<&mut [f64]>,
) -> f64 {
    let mut maxscore = 0.0;
    let ncuts = cands.len();

    // sol != NULL && dircutoffdistweight > 0.0
    if has_incumbent && dircutoffdistweight > 0.0 {
        let mut rng = rng;
        for i in 0..ncuts {
            let cand = &cands[i];
            let mut score;

            let intsupport = if intsupportweight > 0.0 {
                intsupportweight * cand.numintcols / cand.nnz
            } else {
                0.0
            };

            let objparallelism = if objparalweight > 0.0 {
                objparalweight * cand.objpar
            } else {
                0.0
            };

            let efficacy = cand.efficacy;

            // Non-local branch only (pool rows are non-local per)
            let cutoffdist = cand.cutoffdist;
            score = dircutoffdistweight * s4_max(cutoffdist, efficacy);

            let efficacy = efficacy * efficacyweight;
            // C `score += objparallelism + intsupport + efficacy;` (:134): RHS summed first
            score += objparallelism + intsupport + efficacy;

            // add small term to prefer global pool cuts
            if cand.in_global {
                score += 1e-4;
            }

            if let Some(ref mut r) = rng {
                score += r.get_real(0.0, 1e-6);
            }

            maxscore = s4_max(maxscore, score);

            if let Some(ref mut s) = scores {
                s[i] = score;
            }
        }
    } else {
        // in case there is no solution add the directed cutoff distance weight to the efficacy weight
        // since the efficacy underestimates the directed cutoff distance
        let efficacyweight = efficacyweight + dircutoffdistweight;
        let mut rng = rng;
        for i in 0..ncuts {
            let cand = &cands[i];
            let mut score;

            let intsupport = if intsupportweight > 0.0 {
                intsupportweight * cand.numintcols / cand.nnz
            } else {
                0.0
            };

            let objparallelism = if objparalweight > 0.0 {
                objparalweight * cand.objpar
            } else {
                0.0
            };

            let efficacy = if efficacyweight > 0.0 {
                efficacyweight * cand.efficacy
            } else {
                0.0
            };

            score = objparallelism + intsupport + efficacy;

            // add small term to prefer global pool cuts
            if cand.in_global {
                score += 1e-4;
            }

            if let Some(ref mut r) = rng {
                score += r.get_real(0.0, 1e-6);
            }

            maxscore = s4_max(maxscore, score);

            if let Some(ref mut s) = scores {
                s[i] = score;
            }
        }
    }
    maxscore
}

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::selectBestCut (lines 198-227)
pub(crate) fn s4_select_best_cut(perm: &mut [usize], scores: &mut [f64]) {
    let ncuts = perm.len();
    debug_assert!(ncuts > 0);

    let mut bestscore = scores[0];
    let mut bestpos = 0;

    for i in 1..ncuts {
        if scores[i] > bestscore {
            bestpos = i;
            bestscore = scores[i];
        }
    }

    perm.swap(bestpos, 0);
    scores.swap(bestpos, 0);
}

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::filterWithParallelism (lines 231-265)
pub(crate) fn s4_filter_with_parallelism(
    cut: usize,
    perm: &mut [usize],
    scores: &mut [f64],
    ncuts: usize,
    goodscore: f64,
    goodmaxparall: f64,
    maxparall: f64,
    par: &dyn Fn(usize, usize) -> f64,
) -> usize {
    debug_assert!(ncuts == 0 || !perm.is_empty());
    debug_assert!(ncuts == 0 || !scores.is_empty());

    let mut ncuts = ncuts;

    let mut idx = ncuts as i64 - 1;
    while idx >= 0 {
        let i = idx as usize;
        let thisparall = par(cut, perm[i]);
        let thismaxparall = if scores[i] >= goodscore {
            goodmaxparall
        } else {
            maxparall
        };

        if thisparall > thismaxparall {
            ncuts -= 1;
            perm.swap(i, ncuts);
            scores.swap(i, ncuts);
        }

        idx -= 1;
    }

    ncuts
}

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::SCIPselectCutsHybrid (lines 437-522)
pub(crate) fn s4_select_cuts_hybrid(
    all: &[S4Cand],
    nforced: usize,
    rng: Option<&mut S4RandNumGen>,
    has_incumbent: bool,
    goodscorefac: f64,
    badscorefac: f64,
    goodmaxparall: f64,
    maxparall: f64,
    dircutoffdistweight: f64,
    efficacyweight: f64,
    objparalweight: f64,
    intsupportweight: f64,
    maxselected: usize,
    par: &dyn Fn(usize, usize) -> f64,
) -> Vec<usize> {
    debug_assert!(!all.is_empty());
    debug_assert!(nforced <= all.len());

    let mut nselectedcuts = 0usize;

    // (i) score the forced cuts (scores = None; RNG draws included)
    let mut rng = rng;
    let maxforcedscores = s4_scoring(
        &all[..nforced],
        rng.as_deref_mut(),
        dircutoffdistweight,
        efficacyweight,
        objparalweight,
        intsupportweight,
        has_incumbent,
        None,
    );

    // (ii) build perm and scores for non-forced cuts; score them
    let ncuts = all.len() - nforced;
    let mut perm: Vec<usize> = (nforced..all.len()).collect();
    let mut scores: Vec<f64> = vec![0.0; ncuts];
    let maxnonforcedscores = s4_scoring(
        &all[nforced..],
        rng,
        dircutoffdistweight,
        efficacyweight,
        objparalweight,
        intsupportweight,
        has_incumbent,
        Some(&mut scores),
    );

    // (iii) goodscore = max of the two maxima; badscore = goodscore*badfac; goodscore *= goodfac
    let mut goodscore = s4_max(maxforcedscores, maxnonforcedscores);
    let badscore = goodscore * badscorefac;
    goodscore *= goodscorefac;

    // (iv) filter by each forced cut
    let mut ncuts = ncuts;
    for i in 0..nforced {
        if ncuts == 0 {
            break;
        }
        ncuts = s4_filter_with_parallelism(
            i,
            &mut perm,
            &mut scores,
            ncuts,
            goodscore,
            goodmaxparall,
            maxparall,
            par,
        );
    }

    // (v) greedy loop
    let mut base = 0usize;
    while ncuts > 0 {
        s4_select_best_cut(
            &mut perm[base..base + ncuts],
            &mut scores[base..base + ncuts],
        );
        let selectedcut = perm[base];

        // if the best cut of the remaining cuts is considered bad, we discard it and all remaining cuts
        if scores[base] < badscore {
            break;
        }

        nselectedcuts += 1;

        // if the maximal number of cuts was selected, we can stop here
        if nselectedcuts == maxselected {
            break;
        }

        // move the pointers to the next position and filter the remaining cuts
        base += 1;
        ncuts -= 1;

        ncuts = s4_filter_with_parallelism(
            selectedcut,
            &mut perm[base..],
            &mut scores[base..],
            ncuts,
            goodscore,
            goodmaxparall,
            maxparall,
            par,
        );
    }

    perm[..nselectedcuts].to_vec()
}

/// SCIP 11.0.0 (32dd425): cutsel_hybrid.c::cutselSelectHybrid (lines 335-362)
pub(crate) fn s4_cutsel_select_hybrid(
    all: &[S4Cand],
    nforced: usize,
    root: bool,
    rng: Option<&mut S4RandNumGen>,
    has_incumbent: bool,
    maxselected: usize,
    par: &dyn Fn(usize, usize) -> f64,
) -> Vec<usize> {
    let mut minortho = S4_MINORTHO;
    if root {
        minortho = S4_MINORTHOROOT;
    }

    let maxparall = 1.0 - minortho;
    let goodmaxparall = s4_max(0.5, 1.0 - minortho);

    s4_select_cuts_hybrid(
        all,
        nforced,
        rng,
        has_incumbent,
        S4_GOODSCORE,
        S4_BADSCORE,
        goodmaxparall,
        maxparall,
        S4_DIRCUTOFFDISTWEIGHT,
        S4_EFFICACYWEIGHT,
        S4_OBJPARALWEIGHT,
        S4_INTSUPPORTWEIGHT,
        maxselected,
        par,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_cand(
        efficacy: f64,
        objpar: f64,
        intsup_ratio: f64,
        cutoffdist: f64,
        in_global: bool,
    ) -> S4Cand {
        S4Cand {
            pool: 0,
            row: 0,
            efficacy,
            objpar,
            numintcols: intsup_ratio,
            nnz: 1.0,
            cutoffdist,
            in_global,
        }
    }

    #[test]
    fn s4_hybrid_scoring_weight_order() {
        // No RNG, has_incumbent=false, dircutoffdistweight=0.0
        // efficacyweight=1.0, objparalweight=0.1, intsupportweight=0.1
        // else branch: efficacyweight += dircutoffdistweight => 1.0 + 0.0 = 1.0
        // score = objpar*0.1 + intsup*0.1 + eff*1.0
        let cands = vec![
            make_cand(2.0, 3.0, 0.5, 0.0, false),
            make_cand(1.0, 2.0, 0.4, 0.0, true),
        ];
        let mut scores = vec![0.0; 2];
        let maxscore = s4_scoring(
            &cands,
            None,
            0.0,   // dircutoffdistweight
            1.0,   // efficacyweight
            0.1,   // objparalweight
            0.1,   // intsupportweight
            false, // has_incumbent
            Some(&mut scores),
        );

        // cut 0: objpar=0.1*3.0=0.3, intsup=0.1*0.5=0.05, eff=1.0*2.0=2.0
        // score = 0.3 + 0.05 + 2.0 = 2.35
        assert!((scores[0] - 2.35).abs() < 1e-12);

        // cut 1: objpar=0.1*2.0=0.2, intsup=0.1*0.4=0.04, eff=1.0*1.0=1.0
        // score = 0.2 + 0.04 + 1.0 = 1.24, then +1e-4 for in_global => 1.2401
        assert!((scores[1] - 1.2401).abs() < 1e-12);

        assert!((maxscore - 2.35).abs() < 1e-12);
    }

    #[test]
    fn s4_hybrid_incumbent_dircutoff_branch() {
        // has_incumbent=true, dircutoffdistweight=0.5
        // efficacyweight=1.0, objparalweight=0.1, intsupportweight=0.1
        // score = 0.5*max(cutoffdist, eff) + objpar + intsup + eff
        let cands = vec![
            make_cand(2.0, 3.0, 0.5, 1.0, false),
            make_cand(1.0, 2.0, 0.4, 3.0, false),
        ];
        let mut scores = vec![0.0; 2];
        let maxscore = s4_scoring(
            &cands,
            None,
            0.5,  // dircutoffdistweight
            1.0,  // efficacyweight
            0.1,  // objparalweight
            0.1,  // intsupportweight
            true, // has_incumbent
            Some(&mut scores),
        );

        // cut 0: cutoffdist=1.0, eff=2.0 => max=2.0, 0.5*2.0=1.0
        // objpar=0.1*3.0=0.3, intsup=0.1*0.5=0.05, eff=2.0*1.0=2.0
        // score = 1.0 + 0.3 + 0.05 + 2.0 = 3.35
        assert!((scores[0] - 3.35).abs() < 1e-12);

        // cut 1: cutoffdist=3.0, eff=1.0 => max=3.0, 0.5*3.0=1.5
        // objpar=0.1*2.0=0.2, intsup=0.1*0.4=0.04, eff=1.0*1.0=1.0
        // score = 1.5 + 0.2 + 0.04 + 1.0 = 2.74
        assert!((scores[1] - 2.74).abs() < 1e-12);

        assert!((maxscore - 3.35).abs() < 1e-12);
    }

    #[test]
    fn s4_hybrid_incumbent_addition_order_exact() {
        // C :134 `score += objparallelism + intsupport + efficacy` evaluates the RHS first.
        let cands = vec![make_cand(1e-16, 1.0, 0.0, 1.0, false)];
        let mut scores = vec![0.0; 1];
        s4_scoring(&cands, None, 1.0, 1.0, 1e-16, 0.0, true, Some(&mut scores));
        let base = 1.0 * s4_max(1.0, 1e-16);
        let rhs = 1e-16 * 1.0 + 0.0 + 1e-16 * 1.0;
        assert_eq!(scores[0].to_bits(), (base + rhs).to_bits());
        assert_ne!(
            scores[0].to_bits(),
            (((base + 1e-16) + 0.0) + 1e-16).to_bits()
        );
    }

    #[test]
    fn s4_hybrid_scoring_rng_delta() {
        // One get_real(0, 1e-6) draw per cut, in array order.
        let cands = vec![
            make_cand(2.0, 0.0, 0.0, 0.0, false),
            make_cand(1.0, 0.0, 0.0, 0.0, false),
        ];
        let mut scores = vec![0.0; 2];
        let mut rng = S4RandNumGen::new(S4_RANDSEED);
        s4_scoring(
            &cands,
            Some(&mut rng),
            0.0,
            1.0,
            0.0,
            0.0,
            false,
            Some(&mut scores),
        );
        let mut ref_rng = S4RandNumGen::new(S4_RANDSEED);
        let d0 = ref_rng.get_real(0.0, 1e-6);
        let d1 = ref_rng.get_real(0.0, 1e-6);
        assert!((0.0..=1e-6).contains(&d0) && (0.0..=1e-6).contains(&d1));
        assert_eq!(scores[0].to_bits(), (2.0f64 + d0).to_bits());
        assert_eq!(scores[1].to_bits(), (1.0f64 + d1).to_bits());
        assert_eq!(rng.get_rand(), ref_rng.get_rand());
    }

    #[test]
    fn s4_hybrid_else_branch_adds_dircutoff_to_efficacy() {
        // has_incumbent=false, dircutoffdistweight=0.5, efficacyweight=1.0
        // else branch: efficacyweight becomes 1.0 + 0.5 = 1.5
        let cands = vec![make_cand(2.0, 0.0, 0.0, 0.0, false)];
        let mut scores = vec![0.0; 1];
        let maxscore = s4_scoring(
            &cands,
            None,
            0.5,   // dircutoffdistweight
            1.0,   // efficacyweight
            0.0,   // objparalweight
            0.0,   // intsupportweight
            false, // has_incumbent
            Some(&mut scores),
        );

        // efficacyweight = 1.0 + 0.5 = 1.5
        // efficacy = 1.5 * 2.0 = 3.0
        // score = 0.0 + 0.0 + 3.0 = 3.0
        assert!((scores[0] - 3.0).abs() < 1e-12);
        assert!((maxscore - 3.0).abs() < 1e-12);
    }

    #[test]
    fn s4_hybrid_select_best_takes_first_max() {
        // Two equal maxima; the first position wins (strict '>' comparison)
        let mut perm = vec![10, 20, 30, 40];
        let mut scores = vec![5.0, 5.0, 3.0, 1.0];
        s4_select_best_cut(&mut perm, &mut scores);
        // bestpos=0 (first max at index 0), swap 0 with 0 => no change
        assert_eq!(perm[0], 10);
        assert!((scores[0] - 5.0).abs() < 1e-12);

        // Now test where first max is at index 1
        let mut perm2 = vec![10, 20, 30, 40];
        let mut scores2 = vec![3.0, 5.0, 5.0, 1.0];
        s4_select_best_cut(&mut perm2, &mut scores2);
        // bestpos=1 (first max), swap 1 with 0
        assert_eq!(perm2[0], 20);
        assert!((scores2[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn s4_hybrid_filter_from_end_swaps() {
        // cut=0, perm=[1,2,3], scores=[0.5, 0.9, 0.3], ncuts=3
        // goodscore=0.8, goodmaxparall=0.5, maxparall=0.1
        // par(cut, perm[i]): par(0,1)=0.2, par(0,2)=0.6, par(0,3)=0.05
        // Loop i=2: par(0,3)=0.05, scores[2]=0.3 < 0.8 => maxparall=0.1, 0.05 > 0.1? No
        // Loop i=1: par(0,2)=0.6, scores[1]=0.9 >= 0.8 => goodmaxparall=0.5, 0.6 > 0.5? Yes
        //   ncuts=2, swap perm[1] with perm[2], scores[1] with scores[2]
        //   perm=[1,3,2], scores=[0.5,0.3,0.9]
        // Loop i=0: par(0,1)=0.2, scores[0]=0.5 < 0.8 => maxparall=0.1, 0.2 > 0.1? Yes
        //   ncuts=1, swap perm[0] with perm[1], scores[0] with scores[1]
        //   perm=[3,1,2], scores=[0.3,0.5,0.9]
        // Return ncuts=1
        let mut perm = vec![1, 2, 3];
        let mut scores = vec![0.5, 0.9, 0.3];
        let par_fn = |a: usize, b: usize| -> f64 {
            match (a, b) {
                (0, 1) => 0.2,
                (0, 2) => 0.6,
                (0, 3) => 0.05,
                _ => 0.0,
            }
        };
        let ncuts =
            s4_filter_with_parallelism(0, &mut perm, &mut scores, 3, 0.8, 0.5, 0.1, &par_fn);
        assert_eq!(ncuts, 1);
        assert_eq!(perm[0], 3);
        assert!((scores[0] - 0.3).abs() < 1e-12);
    }

    #[test]
    fn s4_hybrid_greedy_stops_at_maxselected() {
        // 3 non-forced cuts, maxselected=2
        // All have high scores, no parallelism filtering
        let cands = vec![
            make_cand(3.0, 0.0, 0.0, 0.0, false),
            make_cand(2.0, 0.0, 0.0, 0.0, false),
            make_cand(1.0, 0.0, 0.0, 0.0, false),
        ];
        let par_fn = |_a: usize, _b: usize| -> f64 { 0.0 };
        let result = s4_select_cuts_hybrid(
            &cands, 0, None, false, 0.9, // goodscorefac
            0.0, // badscorefac
            0.5, // goodmaxparall
            0.1, // maxparall
            0.0, // dircutoffdistweight
            1.0, // efficacyweight
            0.0, // objparalweight
            0.0, // intsupportweight
            2,   // maxselected
            &par_fn,
        );
        assert_eq!(result.len(), 2);
        // First selected should be index 0 (highest efficacy=3.0), second index 1
        assert_eq!(result[0], 0);
        assert_eq!(result[1], 1);
    }

    #[test]
    fn s4_hybrid_badscore_break() {
        // badscorefac=0.5, goodscorefac=1.0
        // Scores: cut0=10, cut1=1, cut2=0.5
        // goodscore = max(0, 10) = 10, badscore = 10*0.5 = 5
        // Greedy: select cut0 (score 10 >= 5), then cut1 (score 1 < 5) => break
        let cands = vec![
            make_cand(10.0, 0.0, 0.0, 0.0, false),
            make_cand(1.0, 0.0, 0.0, 0.0, false),
            make_cand(0.5, 0.0, 0.0, 0.0, false),
        ];
        let par_fn = |_a: usize, _b: usize| -> f64 { 0.0 };
        let result = s4_select_cuts_hybrid(
            &cands, 0, None, false, 1.0, // goodscorefac
            0.5, // badscorefac
            0.5, // goodmaxparall
            0.1, // maxparall
            0.0, // dircutoffdistweight
            1.0, // efficacyweight
            0.0, // objparalweight
            0.0, // intsupportweight
            100, // maxselected
            &par_fn,
        );
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 0);
    }

    #[test]
    fn s4_hybrid_forced_filter() {
        // 1 forced cut (index 0), 2 non-forced cuts (indices 1,2)
        // Forced cut filters out non-forced cut 1 due to parallelism
        // par(0,1)=0.9 > maxparall=0.1 => filtered
        // par(0,2)=0.05 <= maxparall=0.1 => kept
        let cands = vec![
            make_cand(5.0, 0.0, 0.0, 0.0, false), // forced
            make_cand(3.0, 0.0, 0.0, 0.0, false), // non-forced, parallel to forced
            make_cand(2.0, 0.0, 0.0, 0.0, false), // non-forced, not parallel
        ];
        let par_fn = |a: usize, b: usize| -> f64 {
            match (a, b) {
                (0, 1) => 0.9,
                (0, 2) => 0.05,
                _ => 0.0,
            }
        };
        let result = s4_select_cuts_hybrid(
            &cands, 1, // nforced
            None, false, 0.9, // goodscorefac
            0.0, // badscorefac
            0.5, // goodmaxparall
            0.1, // maxparall
            0.0, // dircutoffdistweight
            1.0, // efficacyweight
            0.0, // objparalweight
            0.0, // intsupportweight
            100, // maxselected
            &par_fn,
        );
        // Only cut 2 survives forced filtering
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 2);
    }

    #[test]
    fn s4_hybrid_goodmaxparall_0p5_vs_maxparall_0p1() {
        // root=false: minortho=0.90, maxparall=1.0-0.90=0.1, goodmaxparall=max(0.5,0.1)=0.5
        // goodscorefac=0.9, badscorefac=0.0
        // Two non-forced cuts: cut0 has high score (good), cut1 has low score (not good)
        // par(selected, cut0)=0.3, par(selected, cut1)=0.3
        // For good cut: 0.3 > 0.5? No => kept
        // For non-good cut: 0.3 > 0.1? Yes => removed
        // But we need a selected cut first. Let's use 3 cuts:
        // cut0: efficacy=10 (will be selected first)
        // cut1: efficacy=9 (good score, par=0.3 from cut0 => kept since 0.3 <= 0.5)
        // cut2: efficacy=1 (bad score relative to goodscore, par=0.3 from cut0 => removed since 0.3 > 0.1)
        // goodscore = max(0, 10) = 10, goodscore *= 0.9 => 9.0
        // badscore = 10 * 0.0 = 0.0
        // Greedy: select cut0 (score 10), filter: cut1 score=9 >= 9.0 => goodmaxparall=0.5, par=0.3 <= 0.5 => kept
        //          cut2 score=1 < 9.0 => maxparall=0.1, par=0.3 > 0.1 => removed
        // Then select cut1 (score 9), done (ncuts becomes 0 after filter)
        let cands = vec![
            make_cand(10.0, 0.0, 0.0, 0.0, false),
            make_cand(9.0, 0.0, 0.0, 0.0, false),
            make_cand(1.0, 0.0, 0.0, 0.0, false),
        ];
        let par_fn = |a: usize, b: usize| -> f64 {
            match (a, b) {
                (0, 1) => 0.3,
                (0, 2) => 0.3,
                (1, 2) => 0.0,
                _ => 0.0,
            }
        };
        let result = s4_cutsel_select_hybrid(
            &cands, 0, false, // root=false
            None, false, // has_incumbent
            100,   // maxselected
            &par_fn,
        );
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], 0);
        assert_eq!(result[1], 1);
    }
}
