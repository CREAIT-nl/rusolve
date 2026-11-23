//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! propagateOrbitope, Propagate, AddOrbitope, Reset, GetDefaultColumnOrdering.
//!
//! No gate site exists in this job: nothing else in the crate calls these functions yet,
//! so no `if` gate is needed.

use super::*;
use crate::prelude::*;
use std::sync::Arc;

use super::orbitopal_faces::static_orbitope_fix;
use super::orbitopal_path::{get_column_order, get_row_order};

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_RowOrdering (lines 59-63)
const S3_ROW_ORDERING_NONE: HighsInt = 0;

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateOrbitope (lines 1997-2079)
///
/// Without the debug hash block (lines 2034-2066) and without lastnodenumber (R5).
/// No probing or strong-dual checks (Q6 and Q11).
pub(crate) fn propagate_orbitope(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    orbidata: &S3OrbitopeData,
    orbi_idx: HighsInt,
    path: Option<&Arc<S3PathNode>>,
    is_integral: &[bool],
    eps: f64,
    feastol: f64,
) -> (bool, HighsInt) {
    let mut nfixedvars: HighsInt = 0;
    let mut infeasible = false;

    debug_assert!(orbidata.ncols > 0);
    debug_assert!(orbidata.nrows > 0);

    // get row order
    let roworder: Option<Vec<HighsInt>>;
    let nselrows: HighsInt;

    if orbidata.rowordering == S3_ROW_ORDERING_NONE {
        // C lines 580-585: rowordering NONE early return
        roworder = None;
        nselrows = orbidata.nrows;
    } else {
        let ro = get_row_order(orbi_idx, path);
        nselrows = ro.len() as HighsInt;
        roworder = Some(ro);
    }

    debug_assert!(nselrows >= 0);
    debug_assert!(nselrows <= orbidata.nrows);

    if nselrows == 0 {
        // C lines 2028-2029: goto FREEROWS with nfixedvars=0, infeasible=FALSE
        return (false, 0);
    }

    // get column order
    let colorder = get_column_order(orbi_idx, path, orbidata.columnordering, orbidata.ncols);

    static_orbitope_fix(
        ms,
        localdom,
        orbidata,
        roworder.as_deref(),
        nselrows,
        colorder.as_ref().map(|co| co.0.as_slice()),
        is_integral,
        eps,
        feastol,
    )
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::SCIPorbitopalReductionPropagate (lines 2139-2206)
///
/// Minus the SCIPallowStrongDualReds and SCIPinProbing early returns (Q6 and Q11),
/// minus the debug message.
pub(crate) fn orbitopal_reduction_propagate(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    data: &S3SymData,
    path: Option<&Arc<S3PathNode>>,
    is_integral: &[bool],
    eps: f64,
    feastol: f64,
    stats: &mut S3Stats,
) -> (bool, HighsInt, bool) {
    let mut infeasible = false;
    let mut nred: HighsInt = 0;
    let mut didrun = false;

    // early termination (C lines 2166-2167)
    if data.orbitopes.is_empty() {
        return (false, 0, false);
    }

    // propagate all orbitopes (C lines 2181-2199)
    for orbi_idx in 0..data.orbitopes.len() {
        let (infeas, nfixed) = propagate_orbitope(
            ms,
            localdom,
            &data.orbitopes[orbi_idx],
            orbi_idx as HighsInt,
            path,
            is_integral,
            eps,
            feastol,
        );
        infeasible = infeas;
        nred += nfixed;

        // a symmetry propagator has ran, so set didrun to TRUE (C line 2191)
        didrun = true;

        // stop if we find infeasibility in one of the methods (C lines 2194-2198)
        if infeasible {
            break;
        }
    }

    // C lines 2201-2203
    stats.orbitope_nred += nred as u64;
    if infeasible {
        stats.orbitope_ncutoff += 1;
    }

    (infeasible, nred, didrun)
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::addOrbitope (lines 1097-1247)
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::SCIPorbitopalReductionAddOrbitope (lines 2209-2240)
///
/// Data build only; there are no SCIP var captures, hash maps or events in this port
/// (the path replaces nodeinfos, D3; no nonlinear conshdlr exists in HiGHS, so
/// rowIsBranchRow :982-1017 reduces to the integrality of the row's first variable,
/// via vartypeIsBranchRowType :165-185 with the nonlinear branch always false).
pub(crate) fn orbitopal_reduction_add_orbitope(
    data: &mut S3SymData,
    rowordering: HighsInt,
    colordering: HighsInt,
    vars: &[HighsInt],
    nrows: HighsInt,
    ncols: HighsInt,
    is_integral: &[bool],
) -> bool {
    debug_assert!(nrows >= 0);
    debug_assert!(ncols >= 0);

    let nelem = nrows * ncols;
    debug_assert!(nelem >= 0);

    // prevent trivial case with empty orbitope (C lines 1126-1129)
    if nelem == 0 {
        return false;
    }

    debug_assert!(vars.len() == nelem as usize);

    // build nested vars[row][col] (C lines 1156-1175)
    let mut orbi_vars: Vec<Vec<HighsInt>> = Vec::with_capacity(nrows as usize);
    for row in 0..nrows {
        let mut row_vars: Vec<HighsInt> = Vec::with_capacity(ncols as usize);
        for col in 0..ncols {
            let idx = row * ncols + col;
            let var = vars[idx as usize];
            debug_assert!((var as usize) < is_integral.len());
            row_vars.push(var);
        }
        orbi_vars.push(row_vars);
    }

    // rowindexmap and colindexmap (C lines 1150-1151, 1179-1182)
    // Vec<HighsInt> of length is_integral.len(), initialised to -1 (sentinel: not in the orbitope)
    let mut rowindexmap: Vec<HighsInt> = vec![-1; is_integral.len()];
    let mut colindexmap: Vec<HighsInt> = vec![-1; is_integral.len()];

    for row in 0..nrows {
        for col in 0..ncols {
            let var = orbi_vars[row as usize][col as usize];
            rowindexmap[var as usize] = row;
            colindexmap[var as usize] = col;
        }
    }

    // count number of branchable rows in orbitope (C lines 1187-1194)
    // rowIsBranchRow reduces to is_integral of the row's first variable
    let mut nbranchrows: HighsInt = 0;
    for i in 0..nrows {
        let first_var = orbi_vars[i as usize][0];
        if is_integral[first_var as usize] {
            nbranchrows += 1;
        }
    }

    // Debug-only (#ifndef NDEBUG) lastnodenumber and dbghash lines (C lines 1141-1144) are dropped (R5).
    // nodeinfos and event code (C lines 1199-1213) are dropped (D3: path replaces nodeinfos).

    let orbidata = S3OrbitopeData {
        vars: orbi_vars,
        nrows,
        ncols,
        nbranchrows,
        rowindexmap,
        colindexmap,
        columnordering: colordering,
        rowordering,
    };

    // add orbitope to orbitopal reduction data (C lines 1240-1242)
    data.orbitopes.push(orbidata);

    true
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::SCIPorbitopalReductionReset (lines 2244-2267)
///
/// Clears the dynamic orbitope list. Does not touch static_orbitopes or lexdatas.
pub(crate) fn orbitopal_reduction_reset(data: &mut S3SymData) {
    data.orbitopes.clear();
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::SCIPorbitopalReductionGetDefaultColumnOrdering (lines 2337-2344)
///
/// Returns the default column ordering (DEFAULT_COLUMNORDERING = MEDIAN at line 119, Q5).
pub(crate) fn orbitopal_reduction_default_column_ordering() -> HighsInt {
    S3_COLUMN_ORDERING
}

#[cfg(test)]
mod tests {
    use super::*;

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
        ms.mipdata_.cliquetable = crate::mip::Slot::new(HighsCliqueTable::new(n));
        let implications = HighsImplications::new(&ms);
        ms.mipdata_.implications = crate::mip::Slot::new(implications);
        ms
    }

    fn make_localdom(n: HighsInt, col_lower: &[f64], col_upper: &[f64]) -> HighsDomain {
        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = col_lower.to_vec();
        localdom.col_upper_ = col_upper.to_vec();
        localdom.col_lower_pos_ = vec![-1; n as usize];
        localdom.col_upper_pos_ = vec![-1; n as usize];
        localdom.changedcolsflags_ = vec![0; n as usize];
        localdom
    }

    #[test]
    fn s3_orb_default_ordering_is_median() {
        assert_eq!(
            orbitopal_reduction_default_column_ordering(),
            S3_COLUMN_ORDERING_MEDIAN
        );
    }

    #[test]
    fn s3_orb_add_orbitope_builds_data() {
        let mut data = S3SymData::default();
        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        let result = orbitopal_reduction_add_orbitope(
            &mut data,
            1, // rowordering branching
            4, // colordering median
            &vars,
            3,
            3,
            &is_integral,
        );
        assert!(result);
        assert_eq!(data.orbitopes.len(), 1);
        let orbi = &data.orbitopes[0];
        assert_eq!(orbi.nbranchrows, 3);
        assert_eq!(orbi.vars[1][2], 5);
        assert_eq!(orbi.rowindexmap[5], 1);
        assert_eq!(orbi.colindexmap[5], 2);
        assert_eq!(orbi.columnordering, 4);
        assert_eq!(orbi.rowordering, 1);
    }

    #[test]
    fn s3_orb_add_orbitope_empty_is_rejected() {
        let mut data = S3SymData::default();
        let vars: Vec<HighsInt> = vec![];
        let is_integral = vec![true; 0];
        let result = orbitopal_reduction_add_orbitope(&mut data, 1, 4, &vars, 0, 3, &is_integral);
        assert!(!result);
        assert!(data.orbitopes.is_empty());
    }

    #[test]
    fn s3_orb_reset_clears_orbitopes() {
        let mut data = S3SymData::default();
        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        orbitopal_reduction_add_orbitope(&mut data, 1, 4, &vars, 3, 3, &is_integral);
        assert_eq!(data.orbitopes.len(), 1);
        orbitopal_reduction_reset(&mut data);
        assert!(data.orbitopes.is_empty());
    }

    #[test]
    fn s3_orb_propagate_fixes_3x3_like_static() {
        // 3x3 domain from orbitopal_faces.rs s3_face2_fixes_forced_column_ubs_3x3
        // lower[0]=1; upper[1],[2],[4],[7]=0; all binary, is_integral all true
        let n: HighsInt = 9;
        let mut ms = n_col_solver(n);
        let mut col_lower = vec![0.0f64; 9];
        let mut col_upper = vec![1.0f64; 9];
        col_lower[0] = 1.0;
        col_upper[1] = 0.0;
        col_upper[2] = 0.0;
        col_upper[4] = 0.0;
        col_upper[7] = 0.0;
        let mut localdom = make_localdom(n, &col_lower, &col_upper);

        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        let mut data = S3SymData::default();
        orbitopal_reduction_add_orbitope(&mut data, 0, 0, &vars, 3, 3, &is_integral);

        let orbidata = &data.orbitopes[0];
        let (infeasible, nfixed) = propagate_orbitope(
            &mut ms,
            &mut localdom,
            orbidata,
            0,
            None,
            &is_integral,
            1e-9,
            1e-6,
        );

        assert!(!infeasible);
        assert_eq!(nfixed, 2);
        assert_eq!(localdom.col_upper_[5], 0.0);
        assert_eq!(localdom.col_upper_[8], 0.0);
        assert_eq!(localdom.col_upper_[3], 1.0);
        assert_eq!(localdom.col_upper_[6], 1.0);
        assert_eq!(localdom.col_lower_[0], 1.0);
    }

    #[test]
    fn s3_orb_propagate_identity_colorder_matches_none() {
        // same instance, but colordering 4 (median) and rowordering 0
        let n: HighsInt = 9;
        let mut ms = n_col_solver(n);
        let mut col_lower = vec![0.0f64; 9];
        let mut col_upper = vec![1.0f64; 9];
        col_lower[0] = 1.0;
        col_upper[1] = 0.0;
        col_upper[2] = 0.0;
        col_upper[4] = 0.0;
        col_upper[7] = 0.0;
        let mut localdom = make_localdom(n, &col_lower, &col_upper);

        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        let mut data = S3SymData::default();
        orbitopal_reduction_add_orbitope(&mut data, 0, 4, &vars, 3, 3, &is_integral);

        let orbidata = &data.orbitopes[0];
        let (infeasible, nfixed) = propagate_orbitope(
            &mut ms,
            &mut localdom,
            orbidata,
            0,
            None,
            &is_integral,
            1e-9,
            1e-6,
        );

        assert!(!infeasible);
        assert_eq!(nfixed, 2);
        assert_eq!(localdom.col_upper_[5], 0.0);
        assert_eq!(localdom.col_upper_[8], 0.0);
        assert_eq!(localdom.col_upper_[3], 1.0);
        assert_eq!(localdom.col_upper_[6], 1.0);
        assert_eq!(localdom.col_lower_[0], 1.0);
    }

    #[test]
    fn s3_orb_propagate_branching_rows_none_is_noop() {
        // same instance with rowordering 1 (branching) and path None
        let n: HighsInt = 9;
        let mut ms = n_col_solver(n);
        let mut col_lower = vec![0.0f64; 9];
        let mut col_upper = vec![1.0f64; 9];
        col_lower[0] = 1.0;
        col_upper[1] = 0.0;
        col_upper[2] = 0.0;
        col_upper[4] = 0.0;
        col_upper[7] = 0.0;
        let mut localdom = make_localdom(n, &col_lower, &col_upper);

        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        let mut data = S3SymData::default();
        orbitopal_reduction_add_orbitope(&mut data, 1, 0, &vars, 3, 3, &is_integral);

        let orbidata = &data.orbitopes[0];
        let (infeasible, nfixed) = propagate_orbitope(
            &mut ms,
            &mut localdom,
            orbidata,
            0,
            None,
            &is_integral,
            1e-9,
            1e-6,
        );

        assert!(!infeasible);
        assert_eq!(nfixed, 0);
        // bounds unchanged
        assert_eq!(localdom.col_upper_[5], 1.0);
        assert_eq!(localdom.col_upper_[8], 1.0);
        assert_eq!(localdom.col_upper_[3], 1.0);
        assert_eq!(localdom.col_upper_[6], 1.0);
        assert_eq!(localdom.col_lower_[0], 1.0);
    }

    #[test]
    fn s3_orb_reduction_propagate_counts() {
        // data with one 3x3 orbitope (rowordering 0, colordering 0)
        let n: HighsInt = 9;
        let mut ms = n_col_solver(n);
        let mut col_lower = vec![0.0f64; 9];
        let mut col_upper = vec![1.0f64; 9];
        col_lower[0] = 1.0;
        col_upper[1] = 0.0;
        col_upper[2] = 0.0;
        col_upper[4] = 0.0;
        col_upper[7] = 0.0;
        let mut localdom = make_localdom(n, &col_lower, &col_upper);

        let vars: Vec<HighsInt> = (0..9).collect();
        let is_integral = vec![true; 9];
        let mut data = S3SymData::default();
        orbitopal_reduction_add_orbitope(&mut data, 0, 0, &vars, 3, 3, &is_integral);

        let mut stats = S3Stats::default();
        let (infeasible, nred, didrun) = orbitopal_reduction_propagate(
            &mut ms,
            &mut localdom,
            &data,
            None,
            &is_integral,
            1e-9,
            1e-6,
            &mut stats,
        );

        assert!(!infeasible);
        assert_eq!(nred, 2);
        assert!(didrun);
        assert_eq!(stats.orbitope_nred, 2);
        assert_eq!(stats.orbitope_ncutoff, 0);

        // With an empty data.orbitopes it returns (false, 0, false)
        let empty_data = S3SymData::default();
        let mut stats2 = S3Stats::default();
        let (infeasible2, nred2, didrun2) = orbitopal_reduction_propagate(
            &mut ms,
            &mut localdom,
            &empty_data,
            None,
            &is_integral,
            1e-9,
            1e-6,
            &mut stats2,
        );
        assert!(!infeasible2);
        assert_eq!(nred2, 0);
        assert!(!didrun2);
    }

    #[test]
    fn s3_orb_reduction_propagate_infeasible_cutoff() {
        // a 1x2 orbitope with vars [0,1], lower [0,1], upper [0,1], is_integral [true,true]
        // rowordering 0, colordering 0
        let n: HighsInt = 2;
        let mut ms = n_col_solver(n);
        let col_lower = vec![0.0f64, 1.0];
        let col_upper = vec![0.0f64, 1.0];
        let mut localdom = make_localdom(n, &col_lower, &col_upper);

        let vars: Vec<HighsInt> = vec![0, 1];
        let is_integral = vec![true; 2];
        let mut data = S3SymData::default();
        orbitopal_reduction_add_orbitope(&mut data, 0, 0, &vars, 1, 2, &is_integral);

        let mut stats = S3Stats::default();
        let (infeasible, _nred, didrun) = orbitopal_reduction_propagate(
            &mut ms,
            &mut localdom,
            &data,
            None,
            &is_integral,
            1e-9,
            1e-6,
            &mut stats,
        );

        assert!(infeasible);
        assert!(didrun);
        assert_eq!(stats.orbitope_ncutoff, 1);
    }
}
