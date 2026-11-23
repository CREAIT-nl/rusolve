//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! s3_in_scope, s3_child_path, s3_propagate.
//!
//! No gate site exists in this job: the call sites are in search/part1.rs and
//! search/part2.rs (later jobs), so no `if` gate is needed here.

use super::*;
use crate::prelude::*;
use std::sync::Arc;

use super::debugsol::{
    announce_once, debug_static_orbitopes, debugsol_active, debugsol_objective, debugsol_values,
    dom_alive, kill_scan,
};
use super::lexred::lexicographic_reduction_propagate;
use super::orbitopal::orbitopal_reduction_propagate;
use super::orbitopal_path::orbitopal_branch_info;

/// SCIP 11.0.0 (32dd425): sym_lexorbred.c::symhdlrPropLexOrbRed (lines 1091-1148)
/// SCIP 11.0.0 (32dd425): sym_rowcol.c::symhdlrPropRowCol (lines 2530-2570)
///
/// Returns true only when depth > 0 (SCIPgetDepth > 0), !in_sb, !in_heur (Q6, SCIP
/// skips probing) and symmetry_dynamic is present.
pub(crate) fn s3_in_scope(
    ms: &HighsMipSolver,
    depth: HighsInt,
    in_sb: bool,
    in_heur: bool,
) -> bool {
    if depth <= 0 {
        return false;
    }
    if in_sb {
        return false;
    }
    if in_heur {
        return false;
    }
    ms.mipdata_.symmetry_dynamic.is_some()
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::eventExecNodeBranched (lines 720-962)
///
/// builds the child path node at the branched (parent) node. Reads localdom
/// bounds ONLY (R3): does not mutate localdom and draws no RNG. Returns None when
/// ms.mipdata_.symmetry_dynamic is None. The caller calls this BEFORE the branching
/// change_bound; `depth` is the child depth (root branching = 1, ).
pub(crate) fn s3_child_path(
    ms: &HighsMipSolver,
    localdom: &HighsDomain,
    parent_path: Option<&Arc<S3PathNode>>,
    branch_col: HighsInt,
    depth: HighsInt,
) -> Option<Arc<S3PathNode>> {
    let data = ms.mipdata_.symmetry_dynamic.as_ref()?;
    debug_assert!(depth > 0);
    let orbi = orbitopal_branch_info(
        data,
        &localdom.col_lower_,
        &localdom.col_upper_,
        ms.mipdata_.epsilon,
        parent_path,
        branch_col,
    );
    Some(Arc::new(S3PathNode {
        parent: parent_path.cloned(),
        branch_col,
        depth,
        orbi,
    }))
}

/// SCIP 11.0.0 (32dd425): sym_lexorbred.c::symhdlrPropLexOrbRed (lines 1091-1148)
/// SCIP 11.0.0 (32dd425): sym_rowcol.c::symhdlrPropRowCol (lines 2530-2570)
///
/// HiGHS static orbitopes, then dynamic orbitopal reduction, then LexRed; if
/// anything was reduced, `localdom.propagate` and repeat (up to S3_MAX_PROP_ROUNDS,
/// Q8, approximating SCIP's REDUCEDDOM re-propagation). Returns true when infeasible.
/// SCIPallowStrongDualReds is always true; SCIPinProbing is excluded by the
/// caller scope.
pub(crate) fn s3_propagate(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    path: Option<&Arc<S3PathNode>>,
    depth: HighsInt,
    stats: &mut S3Stats,
) -> bool {
    if !s3_in_scope(ms, depth, false, false) {
        return false;
    }

    let dbg: Option<&'static [f64]> = debugsol_values();

    let kill_watch = match dbg {
        Some(x) if debugsol_active(ms, x) => Some(x),
        _ => None,
    };

    let kill_watch = match kill_watch {
        Some(x) => x,
        None => return s3_propagate_body(ms, localdom, path, depth, stats),
    };

    let obj = debugsol_objective(ms, kill_watch);
    announce_once(obj);

    // SAFETY: ms.mipdata_.symmetry_dynamic is written only in serial setup and
    // never during search (R5). Workers share it read-only.
    let data: *const S3SymData = &**ms.mipdata_.symmetry_dynamic.as_ref().unwrap();
    let data = unsafe { &*data };

    let tol = ms.options_mip_.mip_feasibility_tolerance;
    let static_orbs = debug_static_orbitopes(ms, data);

    let before = dom_alive(
        data,
        &static_orbs,
        kill_watch,
        &localdom.col_lower_,
        &localdom.col_upper_,
        tol,
    );

    let infeasible = s3_propagate_body(ms, localdom, path, depth, stats);

    let after = if infeasible {
        Some(false)
    } else {
        dom_alive(
            data,
            &static_orbs,
            kill_watch,
            &localdom.col_lower_,
            &localdom.col_upper_,
            tol,
        )
    };

    if before == Some(true) && after == Some(false) {
        kill_scan(ms, data, &static_orbs, kill_watch, localdom);
    }

    infeasible
}

fn s3_propagate_body(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    path: Option<&Arc<S3PathNode>>,
    depth: HighsInt,
    stats: &mut S3Stats,
) -> bool {
    if !s3_in_scope(ms, depth, false, false) {
        return false;
    }

    // SAFETY: ms.mipdata_.symmetry_dynamic is written only in serial setup and
    // never during search (R5). Workers share it read-only.
    let data: *const S3SymData = &**ms.mipdata_.symmetry_dynamic.as_ref().unwrap();
    let data = unsafe { &*data };

    let is_integral: Vec<bool> = mip_model!(ms)
        .integrality_
        .iter()
        .map(|&t| t != HighsVarType::kContinuous)
        .collect();

    let eps = ms.mipdata_.epsilon;
    let feastol = ms.options_mip_.mip_feasibility_tolerance;

    for _round in 0..S3_MAX_PROP_ROUNDS {
        // HiGHS static orbitopes (D6/D8 step 1)
        let nfixed_static = ms.with_symmetries_ref(|ms, sym| {
            let mut n: HighsInt = 0;
            for &idx in &data.static_orbitopes {
                n += sym.orbitopes[idx].orbital_fixing(ms, localdom);
                // as HighsSymmetries::propagate_orbitopes (part1.rs:590-595)
                if localdom.infeasible() {
                    break;
                }
            }
            n
        });
        if localdom.infeasible() {
            return true;
        }

        // Dynamic orbitopes
        let (infeasible, nred_orb, _didrun) = orbitopal_reduction_propagate(
            ms,
            localdom,
            data,
            path,
            &is_integral,
            eps,
            feastol,
            stats,
        );
        if infeasible {
            return true;
        }

        // LexRed
        let mut infeasible_lex = false;
        let mut nred_lex: HighsInt = 0;
        let mut didrun = false;
        lexicographic_reduction_propagate(
            ms,
            localdom,
            &data.lexdatas,
            path,
            stats,
            &mut infeasible_lex,
            &mut nred_lex,
            &mut didrun,
        );
        if infeasible_lex {
            return true;
        }

        let nreds = nfixed_static + nred_orb + nred_lex;
        if nreds <= 0 {
            break;
        }
        if localdom.infeasible() {
            return true;
        }
        localdom.propagate(ms);
        if localdom.infeasible() {
            return true;
        }
    }

    localdom.infeasible()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::Slot;
    use crate::mip::symmetry_dynamic::lexred::lexicographic_reduction_add_permutation;

    /// Same layout as the orbitopal.rs test helper: n columns, bounds [0, 1].
    fn make_localdom(n: HighsInt) -> HighsDomain {
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![0.0; n as usize];
        localdom.col_upper_ = vec![1.0; n as usize];
        localdom.col_lower_pos_ = vec![-1; n as usize];
        localdom.col_upper_pos_ = vec![-1; n as usize];
        localdom.changedcolsflags_ = vec![0; n as usize];
        localdom
    }

    fn n_col_solver(n: HighsInt) -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = n;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; n as usize];
        lp.col_lower_ = vec![0.0; n as usize];
        lp.col_upper_ = vec![1.0; n as usize];
        lp.integrality_ = vec![HighsVarType::kInteger; n as usize];
        lp.a_matrix_.num_col_ = n;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; (n + 1) as usize];
        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.mipdata_.cliquetable = Slot::new(HighsCliqueTable::new(n));
        let implications = HighsImplications::new(&ms);
        ms.mipdata_.implications = Slot::new(implications);
        ms
    }

    #[test]
    fn s3_hook_in_scope_requires_depth_sb_heur_and_data() {
        let mut ms = n_col_solver(3);

        // symmetry_dynamic is None by default
        assert!(!s3_in_scope(&ms, 1, false, false));

        ms.mipdata_.symmetry_dynamic = Some(Box::new(S3SymData::default()));

        assert!(s3_in_scope(&ms, 1, false, false));
        assert!(!s3_in_scope(&ms, 0, false, false));
        assert!(!s3_in_scope(&ms, 1, true, false));
        assert!(!s3_in_scope(&ms, 1, false, true));
    }

    #[test]
    fn s3_hook_child_path_none_without_data() {
        let ms = n_col_solver(3);
        let localdom = make_localdom(3);
        let lower_before = localdom.col_lower_.clone();
        let upper_before = localdom.col_upper_.clone();

        let result = s3_child_path(&ms, &localdom, None, 0, 1);
        assert!(result.is_none());
        assert_eq!(localdom.col_lower_, lower_before);
        assert_eq!(localdom.col_upper_, upper_before);
    }

    #[test]
    fn s3_hook_child_path_reads_localdom_only() {
        let mut ms = n_col_solver(12);
        let nrows: HighsInt = 3;
        let ncols: HighsInt = 4;
        let nvars = (nrows * ncols) as usize;

        let mut vars: Vec<Vec<HighsInt>> = Vec::with_capacity(nrows as usize);
        for i in 0..nrows {
            let row: Vec<HighsInt> = (0..ncols).map(|j| i * ncols + j).collect();
            vars.push(row);
        }

        let mut rowindexmap = vec![-1 as HighsInt; nvars];
        let mut colindexmap = vec![-1 as HighsInt; nvars];
        for v in 0..nvars as HighsInt {
            rowindexmap[v as usize] = v / ncols;
            colindexmap[v as usize] = v % ncols;
        }

        let orbidata = S3OrbitopeData {
            vars,
            nrows,
            ncols,
            nbranchrows: 0,
            rowindexmap,
            colindexmap,
            columnordering: S3_COLUMN_ORDERING_MEDIAN,
            rowordering: 1,
        };

        let mut symdata = S3SymData::default();
        symdata.orbitopes.push(orbidata);
        ms.mipdata_.symmetry_dynamic = Some(Box::new(symdata));

        let localdom = make_localdom(12);
        let lower_before = localdom.col_lower_.clone();
        let upper_before = localdom.col_upper_.clone();

        let result = s3_child_path(&ms, &localdom, None, 5, 1);
        assert!(result.is_some());
        let node = result.unwrap();
        assert!(node.parent.is_none());
        assert_eq!(node.branch_col, 5);
        assert_eq!(node.depth, 1);
        assert_eq!(localdom.col_lower_, lower_before);
        assert_eq!(localdom.col_upper_, upper_before);
    }

    #[test]
    fn s3_hook_propagate_depth_zero_is_noop() {
        let mut ms = n_col_solver(3);
        ms.mipdata_.symmetry_dynamic = Some(Box::new(S3SymData::default()));
        let mut localdom = make_localdom(3);
        let lower_before = localdom.col_lower_.clone();
        let upper_before = localdom.col_upper_.clone();
        let mut stats = S3Stats::default();

        let result = s3_propagate(&mut ms, &mut localdom, None, 0, &mut stats);
        assert!(!result);
        assert_eq!(localdom.col_lower_, lower_before);
        assert_eq!(localdom.col_upper_, upper_before);
    }

    #[test]
    fn s3_hook_propagate_depth_one_no_data_is_noop() {
        let mut ms = n_col_solver(3);
        // symmetry_dynamic is None
        let mut localdom = make_localdom(3);
        let lower_before = localdom.col_lower_.clone();
        let upper_before = localdom.col_upper_.clone();
        let mut stats = S3Stats::default();

        let result = s3_propagate(&mut ms, &mut localdom, None, 1, &mut stats);
        assert!(!result);
        assert_eq!(localdom.col_lower_, lower_before);
        assert_eq!(localdom.col_upper_, upper_before);
    }

    #[test]
    fn s3_hook_propagate_depth_one_lexred_reduces() {
        let mut ms = n_col_solver(2);
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.mipdata_.post_solve_stack.initialize_index_maps(0, 2);
        ms.mipdata_.cliquetable = Slot::new(HighsCliqueTable::new(2));
        let imp = crate::mip::implications::HighsImplications::new(&ms);
        *ms.mipdata_.implications = imp;
        ms.mipdata_.domains = vec![Slot::new(HighsDomain::new(&ms))];

        // Swap permutation (0 1), static LexRed (usedynamicorder = false).
        let mut lexdatas: Vec<S3LexData> = Vec::new();
        let mut symvarmap = std::collections::HashMap::new();
        let mut nsymvars: HighsInt = 0;
        let ok = lexicographic_reduction_add_permutation(
            &mut lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &[0, 1],
            &[1, 0],
            false,
        );
        assert!(ok);
        let mut symdata = S3SymData::default();
        symdata.lexdatas = lexdatas;
        ms.mipdata_.symmetry_dynamic = Some(Box::new(symdata));

        // x0 fixed to 0 forces x1 = 0 under x >=_lex swap(x).
        let mut localdom = make_localdom(2);
        localdom.col_upper_[0] = 0.0;
        let mut stats = S3Stats::default();

        let result = s3_propagate(&mut ms, &mut localdom, None, 1, &mut stats);
        assert!(!result);
        assert_eq!(localdom.col_upper_[1], 0.0);
    }
}
