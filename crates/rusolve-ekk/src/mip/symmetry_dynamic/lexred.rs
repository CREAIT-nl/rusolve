//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! LexRed dispatch, Propagate, AddPermutation, statistics.

use super::*;
use crate::mip::symmetry_dynamic::lexred_data::{
    NodeDepthBranchIndex, get_var_order, lexdata_create, lexred_branch_order_from_path,
};
use crate::mip::symmetry_dynamic::lexred_peek::propagate_static_lexred;
use std::collections::HashMap;
use std::sync::Arc;

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateLexredDynamic (lines 1519-1564)
pub(crate) fn propagate_lexred_dynamic(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdata: &S3LexData,
    nodedepthbranchindices: &HashMap<HighsInt, NodeDepthBranchIndex>,
    branchvars: &[HighsInt],
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
) {
    let nvars = lexdata.vars.len() as HighsInt;

    // SCIP line 1549: SCIPallocBufferArray(scip, &varorder, nvars)
    let mut varorder: Vec<HighsInt> = Vec::with_capacity(nvars as usize);

    // SCIP line 1551-1552: getVarOrder
    get_var_order(lexdata, nodedepthbranchindices, branchvars, &mut varorder);

    let nvarorder = varorder.len() as HighsInt;
    debug_assert!(nvarorder >= 0);
    debug_assert!(nvarorder <= nvars);

    // SCIP line 1557: if ( nvarorder > 0 )
    if nvarorder > 0 {
        // SCIP line 1559: propagateStaticLexred(scip, lexdata, varorder, nvarorder, infeasible, nreductions)
        propagate_static_lexred(
            ms,
            localdom,
            lexdata,
            Some(&varorder),
            nvarorder,
            infeasible,
            nreductions,
        );
    }
    // SCIP line 1561: SCIPfreeBufferArray(scip, &varorder) — Rust drop
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateLexredStatic (lines 1569-1592)
pub(crate) fn propagate_lexred_static(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdata: &S3LexData,
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
) {
    // SCIP line 1585: if ( lexdata->nvars == 0 ) return SCIP_OKAY;
    if lexdata.vars.is_empty() {
        return;
    }

    let nvars = lexdata.vars.len() as HighsInt;

    // SCIP line 1589: propagateStaticLexred(scip, lexdata, NULL, lexdata->nvars, infeasible, nreductions)
    propagate_static_lexred(ms, localdom, lexdata, None, nvars, infeasible, nreductions);
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::propagateLexicographicReductionPerm (lines 1597-1630)
pub(crate) fn propagate_lexicographic_reduction_perm(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdata: &S3LexData,
    nodedepthbranchindices: &HashMap<HighsInt, NodeDepthBranchIndex>,
    branchvars: &[HighsInt],
    infeasible: &mut bool,
    nreductions: &mut HighsInt,
) {
    // SCIP line 1619: *nreductions = 0;
    *nreductions = 0;
    // SCIP line 1620: *infeasible = FALSE;
    *infeasible = false;

    // SCIP line 1622: if ( lexdata->isdynamic )
    if lexdata.isdynamic {
        // SCIP line 1624-1625: propagateLexredDynamic
        propagate_lexred_dynamic(
            ms,
            localdom,
            lexdata,
            nodedepthbranchindices,
            branchvars,
            infeasible,
            nreductions,
        );
    } else {
        // SCIP line 1629: propagateLexredStatic
        propagate_lexred_static(ms, localdom, lexdata, infeasible, nreductions);
    }
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::SCIPlexicographicReductionPropagate (lines 1905-2013)
pub(crate) fn lexicographic_reduction_propagate(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    lexdatas: &[S3LexData],
    path: Option<&Arc<S3PathNode>>,
    stats: &mut S3Stats,
    infeasible: &mut bool,
    nred: &mut HighsInt,
    didrun: &mut bool,
) {
    // SCIP line 1932: *infeasible = FALSE;
    *infeasible = false;
    // SCIP line 1933: *nred = 0;
    *nred = 0;
    // SCIP line 1934: *didrun = FALSE;
    *didrun = false;

    // SCIP line 1937: if ( masterdata->nlexdatas == 0 ) return SCIP_OKAY;
    if lexdatas.is_empty() {
        return;
    }

    // SCIP line 1943: if ( masterdata->hasdynamicperm )
    let hasdynamicperm = lexdatas.iter().any(|l| l.isdynamic);

    let mut nodedepthbranchindices: HashMap<HighsInt, NodeDepthBranchIndex> = HashMap::new();
    let mut branchvars: Vec<HighsInt> = Vec::new();

    if hasdynamicperm {
        // SCIP lines 1957-1960: fill nodedepthbranchindices and branchvars from shadowtree.
        // In HiGHS the path always resolves (no inforesolved check).
        let (ndbi, bvars) = lexred_branch_order_from_path(path);
        nodedepthbranchindices = ndbi;
        // SCIP lines 1712-1716: branchvars holds only variables in masterdata->symvarmap
        // ("ignore variables that are irrelevant for lexicographic reduction"). The global
        // symvarmap keys are exactly the vars of the dynamic lexdatas (lexdataCreate :371-402),
        // i.e. the keys of their varmaps. Keep leaf-to-root first-seen order.
        branchvars = bvars
            .into_iter()
            .filter(|col| {
                lexdatas
                    .iter()
                    .any(|l| l.isdynamic && l.varmap.contains_key(col))
            })
            .collect();
    }

    // SCIP line 1973: if ( nbranchvars > 0 || ! masterdata->hasdynamicperm )
    if !branchvars.is_empty() || !hasdynamicperm {
        // SCIP line 1976: for (p = 0; p < masterdata->nlexdatas; ++p)
        for lexdata in lexdatas.iter() {
            // SCIP line 1980-1981: propagateLexicographicReductionPerm
            let mut nlocalred: HighsInt = 0;
            propagate_lexicographic_reduction_perm(
                ms,
                localdom,
                lexdata,
                &nodedepthbranchindices,
                &branchvars,
                infeasible,
                &mut nlocalred,
            );

            // SCIP line 1984: *nred += nlocalred;
            *nred += nlocalred;

            // SCIP line 1987: *didrun = TRUE;
            *didrun = true;

            // SCIP line 1990: if ( *infeasible ) break;
            if *infeasible {
                break;
            }
        }

        // SCIP line 1995: masterdata->nred += *nred;
        stats.lexred_nred += *nred as u64;
        // SCIP line 1996-1997: if ( *infeasible ) ++masterdata->ncutoff;
        if *infeasible {
            stats.lexred_ncutoff += 1;
        }
    }

    // SCIP lines 2001-2010: undo/clean nodedepthbranchindices — Rust drop handles this.
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::SCIPlexicographicReductionAddPermutation (lines 2017-2077)
///
/// Signed permutations are not supported: no symtype argument.
pub(crate) fn lexicographic_reduction_add_permutation(
    lexdatas: &mut Vec<S3LexData>,
    symvarmap: &mut HashMap<HighsInt, HighsInt>,
    nsymvars: &mut HighsInt,
    permvars: &[HighsInt],
    perm: &[HighsInt],
    usedynamicorder: bool,
) -> bool {
    // SCIP line 2069-2070: lexdataCreate
    let ld = lexdata_create(symvarmap, nsymvars, permvars, perm, usedynamicorder);

    // SCIP line 2073-2074: if ( *success ) ++masterdata->nlexdatas;
    match ld {
        Some(ld) => {
            lexdatas.push(ld);
            true
        }
        None => false,
    }
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::SCIPlexicographicReductionReset (lines 2081-2115)
pub(crate) fn lexicographic_reduction_reset(
    lexdatas: &mut Vec<S3LexData>,
    symvarmap: &mut HashMap<HighsInt, HighsInt>,
    nsymvars: &mut HighsInt,
) {
    // SCIP line 2093-2096: while ( masterdata->nlexdatas > 0 ) lexdataFree(...)
    lexdatas.clear();

    // SCIP line 2099-2101: free lexdas array, set to NULL, maxnlexdatas = 0
    // (Rust: clear already handles the Vec; no explicit capacity reset needed)

    // SCIP line 2105-2110: if ( symvarmap != NULL ) { free; nsymvars = 0; }
    symvarmap.clear();
    *nsymvars = 0;

    // SCIP line 2112: masterdata->hasdynamicperm = FALSE;
    // (hasdynamicperm is derived from lexdatas in propagate; no field to reset)
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::SCIPlexicographicReductionGetStatistics (lines 1854-1869)
pub(crate) fn lexicographic_reduction_get_statistics(stats: &S3Stats) -> (u64, u64) {
    // SCIP line 1865: *nred = masterdata->nred;
    // SCIP line 1866: *ncutoff = masterdata->ncutoff;
    (stats.lexred_nred, stats.lexred_ncutoff)
}

/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::SCIPlexicographicReductionPrintStatistics (lines 1872-1901)
///
/// Returns the text; no printing side effect.
pub(crate) fn lexicographic_reduction_print_statistics(lexdatas: &[S3LexData]) -> String {
    // SCIP line 1882: if ( masterdata->nlexdatas == 0 )
    if lexdatas.is_empty() {
        // SCIP line 1884: "   lexicographic reduction:   no permutations\n"
        return "   lexicographic reduction:   no permutations\n".to_string();
    }

    // SCIP line 1888-1889: "   lexicographic reduction: %4d permutations with support sizes "
    let mut result = format!(
        "   lexicographic reduction: {:4} permutations with support sizes ",
        lexdatas.len()
    );

    // SCIP line 1891-1896: for (i = 0; i < masterdata->nlexdatas; ++i)
    for (i, ld) in lexdatas.iter().enumerate() {
        // SCIP line 1893-1894: if ( i > 0 ) ", "
        if i > 0 {
            result.push_str(", ");
        }
        // SCIP line 1895: "%d", masterdata->lexdatas[i]->nvars
        result.push_str(&ld.vars.len().to_string());
    }

    // SCIP line 1898: "\n"
    result.push('\n');

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helpers copied from lexred_peek.rs tests.
    fn int_solver(num_col: HighsInt) -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = num_col;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; num_col as usize];
        lp.col_lower_ = vec![0.0; num_col as usize];
        lp.col_upper_ = vec![10.0; num_col as usize];
        lp.integrality_ = vec![HighsVarType::kInteger; num_col as usize];
        lp.a_matrix_.num_col_ = num_col;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; num_col as usize + 1];
        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0)
    }

    fn make_domain(lb0: f64, ub0: f64, lb1: f64, ub1: f64) -> HighsDomain {
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![lb0, lb1];
        localdom.col_upper_ = vec![ub0, ub1];
        localdom.col_lower_pos_ = vec![-1; 2];
        localdom.col_upper_pos_ = vec![-1; 2];
        localdom.changedcolsflags_ = vec![0; 2];
        localdom
    }

    // Attachments from s3_peek_static_propagation_tightens_lb: a tightening of an
    // integer column runs add_implications, so the clique table and domain are needed.
    fn attach_domain(ms: &mut HighsMipSolver) {
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(0, 2);
        ms.mipdata_.cliquetable = crate::mip::Slot::new(HighsCliqueTable::new(2));
        let imp = crate::mip::implications::HighsImplications::new(ms);
        *ms.mipdata_.implications = imp;
        ms.mipdata_.domains = vec![crate::mip::Slot::new(HighsDomain::new(ms))];
    }

    fn swap_perm(usedynamicorder: bool) -> (Vec<S3LexData>, HashMap<HighsInt, HighsInt>, HighsInt) {
        let mut lexdatas = Vec::new();
        let mut symvarmap = HashMap::new();
        let mut nsymvars: HighsInt = 0;
        let ok = lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1],
            &[1, 0],
            usedynamicorder,
        );
        assert!(ok);
        (lexdatas, symvarmap, nsymvars)
    }

    fn branch_path(col: HighsInt) -> Arc<S3PathNode> {
        Arc::new(S3PathNode {
            parent: None,
            branch_col: col,
            depth: 1,
            orbi: vec![],
        })
    }

    /// Classic LexRed case: swap of two columns, branch x1 = 0 (dynamic order, rooted path)
    /// gives x2 ub = 0.
    #[test]
    fn s3_lxr_swap_branch_tightens_partner_ub() {
        let mut ms = int_solver(2);
        attach_domain(&mut ms);
        let mut localdom = make_domain(0.0, 1.0, 0.0, 1.0);
        let (lexdatas, _symvarmap, _nsymvars) = swap_perm(true);
        localdom.col_upper_[0] = 0.0;
        let path = branch_path(0);
        let mut stats = S3Stats::default();
        let mut infeasible = false;
        let mut nred: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            &mut ms,
            &mut localdom,
            &lexdatas,
            Some(&path),
            &mut stats,
            &mut infeasible,
            &mut nred,
            &mut didrun,
        );
        assert!(!infeasible);
        assert!(didrun);
        assert_eq!(nred, 1);
        assert_eq!(localdom.col_upper_[1], 0.0);
        assert_eq!(stats.lexred_nred, 1);
        assert_eq!(stats.lexred_ncutoff, 0);
    }

    /// Dynamic permutation with no branching on its support: nothing runs.
    #[test]
    fn s3_lxr_dynamic_without_branch_is_noop() {
        let mut ms = int_solver(2);
        attach_domain(&mut ms);
        let mut localdom = make_domain(0.0, 1.0, 0.0, 1.0);
        let (lexdatas, _symvarmap, _nsymvars) = swap_perm(true);
        let mut stats = S3Stats::default();
        let mut infeasible = false;
        let mut nred: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            &mut ms,
            &mut localdom,
            &lexdatas,
            None,
            &mut stats,
            &mut infeasible,
            &mut nred,
            &mut didrun,
        );
        assert!(!didrun);
        assert_eq!(nred, 0);
        assert_eq!(localdom.col_upper_[1], 1.0);
    }

    /// Static permutation propagates without any path.
    #[test]
    fn s3_lxr_static_perm_propagates_without_path() {
        let mut ms = int_solver(2);
        attach_domain(&mut ms);
        let mut localdom = make_domain(0.0, 1.0, 0.0, 1.0);
        let (lexdatas, _symvarmap, _nsymvars) = swap_perm(false);
        localdom.col_upper_[0] = 0.0;
        let mut stats = S3Stats::default();
        let mut infeasible = false;
        let mut nred: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            &mut ms,
            &mut localdom,
            &lexdatas,
            None,
            &mut stats,
            &mut infeasible,
            &mut nred,
            &mut didrun,
        );
        assert!(didrun);
        assert_eq!(localdom.col_upper_[1], 0.0);
    }

    /// Mixed dynamic + static permutations, branching only on a column outside the
    /// symvarmap (support of the static perm): SCIP's nbranchvars is 0, so nothing runs,
    /// not even the static permutation (symmetry_lexred.c:1712-1716, :1973).
    #[test]
    fn s3_lxr_branch_outside_symvars_skips_static() {
        let mut ms = int_solver(4);
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(0, 4);
        ms.mipdata_.cliquetable = crate::mip::Slot::new(HighsCliqueTable::new(4));
        let imp = crate::mip::implications::HighsImplications::new(&ms);
        *ms.mipdata_.implications = imp;
        ms.mipdata_.domains = vec![crate::mip::Slot::new(HighsDomain::new(&ms))];
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![0.0; 4];
        localdom.col_upper_ = vec![1.0; 4];
        localdom.col_lower_pos_ = vec![-1; 4];
        localdom.col_upper_pos_ = vec![-1; 4];
        localdom.changedcolsflags_ = vec![0; 4];
        let mut lexdatas = Vec::new();
        let mut symvarmap = HashMap::new();
        let mut nsymvars: HighsInt = 0;
        assert!(lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1, 2, 3],
            &[1, 0, 2, 3],
            true,
        ));
        assert!(lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1, 2, 3],
            &[0, 1, 3, 2],
            false,
        ));
        assert_eq!(nsymvars, 2);
        localdom.col_upper_[2] = 0.0;
        let path = branch_path(2);
        let mut stats = S3Stats::default();
        let mut infeasible = false;
        let mut nred: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            &mut ms,
            &mut localdom,
            &lexdatas,
            Some(&path),
            &mut stats,
            &mut infeasible,
            &mut nred,
            &mut didrun,
        );
        assert!(!didrun);
        assert!(!infeasible);
        assert_eq!(nred, 0);
        assert_eq!(localdom.col_upper_[3], 1.0);
        assert_eq!(stats.lexred_nred, 0);
    }

    /// Branching x1 = 0 while x2 is forced to 1 is infeasible and counts one cutoff.
    #[test]
    fn s3_lxr_infeasible_counts_cutoff() {
        let mut ms = int_solver(2);
        attach_domain(&mut ms);
        let mut localdom = make_domain(0.0, 1.0, 1.0, 1.0);
        let (lexdatas, _symvarmap, _nsymvars) = swap_perm(true);
        localdom.col_upper_[0] = 0.0;
        let path = branch_path(0);
        let mut stats = S3Stats::default();
        let mut infeasible = false;
        let mut nred: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            &mut ms,
            &mut localdom,
            &lexdatas,
            Some(&path),
            &mut stats,
            &mut infeasible,
            &mut nred,
            &mut didrun,
        );
        assert!(infeasible);
        assert_eq!(stats.lexred_ncutoff, 1);
    }
}
