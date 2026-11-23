//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! lexmin and lexmax face computation for static orbitope propagation.
//! comparison of the two faces and the resulting bound fixings.
//!
//! No gate site exists in this job: nothing else in the crate calls these functions yet,
//! so no `if` gate is needed.

use super::*;
use crate::prelude::*;

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::getArrayEntryOrIndex (lines 512-524)
fn get_array_entry_or_index(arr: Option<&[HighsInt]>, idx: HighsInt) -> HighsInt {
    debug_assert!(idx >= 0);
    if arr.is_none() {
        return idx;
    }
    arr.unwrap()[idx as usize]
}

/// Lexmin and lexmax faces of a static orbitope matrix.
///
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateStaticOrbitope (lines 1479-1869), face computation only
#[derive(Debug, Clone, Default)]
pub(crate) struct S3OrbitopeFaces {
    pub(crate) lexminface: Vec<f64>,
    pub(crate) lexminepsrow: Vec<HighsInt>,
    pub(crate) lexmaxface: Vec<f64>,
    pub(crate) lexmaxepsrow: Vec<HighsInt>,
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateStaticOrbitope (lines 1479-1869), face computation only
pub(crate) fn static_orbitope_faces(
    orbidata: &S3OrbitopeData,
    roworder: Option<&[HighsInt]>,
    nselrows: HighsInt,
    colorder: Option<&[HighsInt]>,
    col_lower: &[f64],
    col_upper: &[f64],
    is_integral: &[bool],
    eps: f64,
) -> Option<S3OrbitopeFaces> {
    debug_assert!(nselrows >= 0);
    debug_assert!(orbidata.nrows > 0);
    debug_assert!(orbidata.nrows >= nselrows);
    let ncols = orbidata.ncols;
    debug_assert!(ncols > 1);

    // nothing to propagate without any rows
    if nselrows <= 0 {
        return Some(S3OrbitopeFaces::default());
    }

    let nelem = nselrows * ncols;

    // compute lexmin face
    let (lexminface, lexminepsrow) = compute_lexmin_face(
        orbidata,
        roworder,
        nselrows,
        colorder,
        col_lower,
        col_upper,
        is_integral,
        eps,
        nelem,
        ncols,
    )?;

    // #ifndef NDEBUG: assertIsOrbitopeMatrix (lines 1712-1715) dropped.

    // compute lexmax face
    let (lexmaxface, lexmaxepsrow) = compute_lexmax_face(
        orbidata,
        roworder,
        nselrows,
        colorder,
        col_lower,
        col_upper,
        is_integral,
        eps,
        nelem,
        ncols,
    )?;

    // #ifndef NDEBUG: assertIsOrbitopeMatrix (lines 1871-1874) dropped.

    Some(S3OrbitopeFaces {
        lexminface,
        lexminepsrow,
        lexmaxface,
        lexmaxepsrow,
    })
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateStaticOrbitope (lines 1560-1710), lexmin face
fn compute_lexmin_face(
    orbidata: &S3OrbitopeData,
    roworder: Option<&[HighsInt]>,
    nselrows: HighsInt,
    colorder: Option<&[HighsInt]>,
    col_lower: &[f64],
    col_upper: &[f64],
    is_integral: &[bool],
    eps: f64,
    nelem: HighsInt,
    ncols: HighsInt,
) -> Option<(Vec<f64>, Vec<HighsInt>)> {
    let mut lexminface = vec![0.0f64; nelem as usize];

    // array to store for each column at which row we add an infinitesimal value, initially at none (-1)
    let mut lexminepsrow = vec![-1 as HighsInt; ncols as usize];
    for colid in 0..ncols {
        lexminepsrow[colid as usize] = -1;
    }

    // last column takes the minimally possible values.
    let colid = ncols - 1;
    let origcolid = get_array_entry_or_index(colorder, colid);
    let mut rowid: HighsInt = 0;
    let mut i = colid;
    while rowid < nselrows {
        let origrowid = get_array_entry_or_index(roworder, rowid);
        let origidx = origrowid * ncols + origcolid;
        let var = orbidata.vars[origrowid as usize][origcolid as usize];

        debug_assert!(i == rowid * ncols + colid);
        debug_assert!(var >= 0);

        lexminface[i as usize] = col_lower[var as usize];

        rowid += 1;
        i += ncols;
    }

    // all previous columns: One-column replacement algorithm
    let mut colid = ncols - 2;
    while colid >= 0 {
        // "rowid" of the last unfixed entry whose domain allows for larger values than the current chosen.
        // If there is none, -1.
        let mut lastunfixed: HighsInt = -1;
        // whether column "colid" is the same as column "colid + 1" up (but excluding) to "rowid"
        let mut iseq = true;

        let origcolid = get_array_entry_or_index(colorder, colid);
        let mut rowid: HighsInt = 0;
        let mut i = colid;
        while rowid < nselrows {
            let origrowid = get_array_entry_or_index(roworder, rowid);
            let origidx = origrowid * ncols + origcolid;
            debug_assert!(i == rowid * ncols + colid);

            // the entry one to the right is not the first column
            debug_assert!((i + 1) % ncols > 0);

            let var = orbidata.vars[origrowid as usize][origcolid as usize];
            debug_assert!(var >= 0);

            if iseq {
                // equality holds up to this row
                // Compare to the entry value on the column immediately right.
                // The value we choose on the left must be at least this.
                // 2 Options:
                // Option 1: The upper bound is smaller. Then we're in an infeasible situation. Resolve as described below.
                // Option 2: The upper bound is greater or equal.
                let ub = col_upper[var as usize];

                // compare to the value in the column right of it
                if sym_lt(ub, lexminface[(i + 1) as usize], eps, kHighsInf)
                    || (lexminepsrow[(colid + 1) as usize] == rowid
                        && sym_eq(ub, lexminface[(i + 1) as usize], eps, kHighsInf))
                {
                    // value of this column can only be strictly smaller than the value in the column to its right
                    // This may not be possible.
                    // Try to repair: Go back to the last row with "room" left, and make the value minimally larger.
                    if lastunfixed >= 0 {
                        // repair: return to the last row with "room", and increase the lexmin-value at that row.
                        debug_assert!(sym_eq(
                            lexminface[(lastunfixed * ncols + colid) as usize],
                            lexminface[(lastunfixed * ncols + colid + 1) as usize],
                            eps,
                            kHighsInf
                        ));
                        let othervar = orbidata.vars
                            [get_array_entry_or_index(roworder, lastunfixed) as usize]
                            [origcolid as usize];
                        if is_integral[othervar as usize] {
                            // discrete type with unit steps: Add one to the bound.
                            debug_assert!(
                                lexminface[(lastunfixed * ncols + colid) as usize]
                                    == (lexminface[(lastunfixed * ncols + colid) as usize] as i64)
                                        as f64
                            );
                            lexminface[(lastunfixed * ncols + colid) as usize] += 1.0;
                            debug_assert!(
                                lexminface[(lastunfixed * ncols + colid) as usize]
                                    == (lexminface[(lastunfixed * ncols + colid) as usize] as i64)
                                        as f64
                            );
                            debug_assert!(sym_le(
                                lexminface[(lastunfixed * ncols + colid) as usize],
                                col_upper[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                        } else {
                            // continuous type, so add an infinitesimal value to the bound
                            debug_assert!(sym_le(
                                lexminface[(lastunfixed * ncols + colid) as usize],
                                col_upper[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                            debug_assert!(lexminepsrow[colid as usize] == -1);
                            lexminepsrow[colid as usize] = lastunfixed;
                        }
                        // now row "lastunfixed" is greater. Restart from here.
                        iseq = false;
                        // the next iteration considers "lastunfixed + 1"
                        rowid = lastunfixed + 1;
                        i = rowid * ncols + colid;
                        continue;
                    } else {
                        // cannot repair. It is infeasible.
                        // SCIPdebugMessage dropped.
                        return None;
                    }
                } else {
                    debug_assert!(sym_ge(ub, lexminface[(i + 1) as usize], eps, kHighsInf));
                    let lb = col_lower[var as usize];
                    debug_assert!(sym_le(lb, ub, eps, kHighsInf));
                    lexminface[i as usize] = cpp_max(lexminface[(i + 1) as usize], lb);
                    debug_assert!(sym_ge(
                        lexminface[i as usize],
                        lexminface[(i + 1) as usize],
                        eps,
                        kHighsInf
                    ));

                    // are we still equal?
                    if sym_gt(
                        lexminface[i as usize],
                        lexminface[(i + 1) as usize],
                        eps,
                        kHighsInf,
                    ) {
                        iseq = false;
                    } else if lexminepsrow[(colid + 1) as usize] == rowid {
                        debug_assert!(sym_eq(
                            lexminface[i as usize],
                            lexminface[(i + 1) as usize],
                            eps,
                            kHighsInf
                        ));
                        debug_assert!(
                            !is_integral
                                [orbidata.vars[origrowid as usize][origcolid as usize] as usize]
                        );
                        debug_assert!(!is_integral[var as usize]);
                        // right column (colid+1) has value x + epsilon, left column (colid) has value x, now
                        // must also become  x + epsilon in order to be larger or equal
                        // by axioms, we can squeeze infinitesimals between one other; epsilon > epsilon.
                        iseq = false;
                        debug_assert!(lexminepsrow[colid as usize] == -1);
                        lexminepsrow[colid as usize] = rowid;
                    }

                    // is there room left to increase this variable further?
                    if is_integral[var as usize] {
                        // discrete type with unit steps: Add one to the bound.
                        debug_assert!(
                            lexminface[i as usize] == (lexminface[i as usize] as i64) as f64
                        );
                        if sym_le(lexminface[i as usize] + 1.0, ub, eps, kHighsInf) {
                            lastunfixed = rowid;
                        }
                    } else {
                        // continuous type: if we can add an infinitesimal value to the current lexminface[i] value,
                        // mark row as 'lastunfixed'
                        if sym_lt(lexminface[i as usize], ub, eps, kHighsInf) {
                            lastunfixed = rowid;
                        }
                    }
                }
            } else {
                // there had been a row before which breaks the equality-condition, choose minimally possible value
                lexminface[i as usize] = col_lower[var as usize];
            }

            rowid += 1;
            i += ncols;
        }

        colid -= 1;
    }

    Some((lexminface, lexminepsrow))
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateStaticOrbitope (lines 1717-1869), lexmax face
fn compute_lexmax_face(
    orbidata: &S3OrbitopeData,
    roworder: Option<&[HighsInt]>,
    nselrows: HighsInt,
    colorder: Option<&[HighsInt]>,
    col_lower: &[f64],
    col_upper: &[f64],
    is_integral: &[bool],
    eps: f64,
    nelem: HighsInt,
    ncols: HighsInt,
) -> Option<(Vec<f64>, Vec<HighsInt>)> {
    let mut lexmaxface = vec![0.0f64; nelem as usize];

    // array to store for each column at which row we subtract an infinitesimal value, initially at none (-1)
    let mut lexmaxepsrow = vec![-1 as HighsInt; ncols as usize];
    for colid in 0..ncols {
        lexmaxepsrow[colid as usize] = -1;
    }

    // first column, fill all unfixed entries with maximally possible values
    let colid = 0;
    let origcolid = get_array_entry_or_index(colorder, colid);
    let mut rowid: HighsInt = 0;
    let mut i = colid;
    while rowid < nselrows {
        let origrowid = get_array_entry_or_index(roworder, rowid);
        let origidx = origrowid * ncols + origcolid;
        let var = orbidata.vars[origrowid as usize][origcolid as usize];

        debug_assert!(i == rowid * ncols + colid);
        debug_assert!(var >= 0);

        lexmaxface[i as usize] = col_upper[var as usize];

        rowid += 1;
        i += ncols;
    }

    // all next columns: One-column replacement algorithm
    let mut colid: HighsInt = 1;
    while colid < ncols {
        // "rowid" of the last unfixed entry whose domain allows for smaller values than the current chosen.
        // If there is none, -1.
        let mut lastunfixed: HighsInt = -1;
        // whether column "colid" is the same as column "colid - 1" up (but excluding) to "rowid"
        let mut iseq = true;

        let origcolid = get_array_entry_or_index(colorder, colid);
        let mut rowid: HighsInt = 0;
        let mut i = colid;
        while rowid < nselrows {
            let origrowid = get_array_entry_or_index(roworder, rowid);
            let origidx = origrowid * ncols + origcolid;
            debug_assert!(i == rowid * ncols + colid);

            // the entry one to the left is not the last column, i.e., this column cannot be the first column
            debug_assert!(i % ncols > 0);

            let var = orbidata.vars[origrowid as usize][origcolid as usize];
            debug_assert!(var >= 0);

            if iseq {
                // equality holds up to this row
                // Compare to the entry value on the column immediately left.
                // The value we choose on the right must be at most this.
                // 2 Options:
                // Option 1: The lower bound is larger. Then we're in an infeasible situation. Resolve as described below.
                // Option 2: The lower bound is smaller or equal.
                let lb = col_lower[var as usize];

                // compare to the value in the column left of it
                if sym_gt(lb, lexmaxface[(i - 1) as usize], eps, kHighsInf)
                    || (lexmaxepsrow[(colid - 1) as usize] == rowid
                        && sym_eq(lb, lexmaxface[(i - 1) as usize], eps, kHighsInf))
                {
                    // value of this column can only be strictly larger than the value in the column to its left
                    // This may not be possible.
                    // Try to repair: Go back to the last row with "room" left, and make the value minimally smaller.
                    if lastunfixed >= 0 {
                        // repair: return to the last row with "room", and decrease the lexmax-value at that row.
                        debug_assert!(sym_eq(
                            lexmaxface[(lastunfixed * ncols + colid) as usize],
                            lexmaxface[(lastunfixed * ncols + colid - 1) as usize],
                            eps,
                            kHighsInf
                        ));
                        let othervar = orbidata.vars
                            [get_array_entry_or_index(roworder, lastunfixed) as usize]
                            [origcolid as usize];
                        if is_integral[othervar as usize] {
                            // discrete type with unit steps: Remove one from the lexmax-value.
                            debug_assert!(
                                lexmaxface[(lastunfixed * ncols + colid) as usize]
                                    == (lexmaxface[(lastunfixed * ncols + colid) as usize] as i64)
                                        as f64
                            );
                            lexmaxface[(lastunfixed * ncols + colid) as usize] -= 1.0;
                            debug_assert!(
                                lexmaxface[(lastunfixed * ncols + colid) as usize]
                                    == (lexmaxface[(lastunfixed * ncols + colid) as usize] as i64)
                                        as f64
                            );
                            debug_assert!(sym_ge(
                                lexmaxface[(lastunfixed * ncols + colid) as usize],
                                col_lower[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                            debug_assert!(sym_le(
                                lexmaxface[(lastunfixed * ncols + colid) as usize],
                                col_upper[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                        } else {
                            // continuous type, so subtract an infinitesimal value to the bound
                            debug_assert!(sym_ge(
                                lexmaxface[(lastunfixed * ncols + colid) as usize],
                                col_lower[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                            debug_assert!(sym_le(
                                lexmaxface[(lastunfixed * ncols + colid) as usize],
                                col_upper[othervar as usize],
                                eps,
                                kHighsInf
                            ));
                            debug_assert!(lexmaxepsrow[colid as usize] == -1);
                            lexmaxepsrow[colid as usize] = lastunfixed;
                        }
                        // now row "lastunfixed" is greater. Restart from here.
                        iseq = false;
                        // the next iteration considers "lastunfixed + 1"
                        rowid = lastunfixed + 1;
                        i = rowid * ncols + colid;
                        continue;
                    } else {
                        // cannot repair. It is infeasible.
                        // SCIPdebugMessage dropped.
                        return None;
                    }
                } else {
                    debug_assert!(sym_le(lb, lexmaxface[(i - 1) as usize], eps, kHighsInf));
                    let ub = col_upper[var as usize];
                    debug_assert!(sym_le(lb, ub, eps, kHighsInf));
                    lexmaxface[i as usize] = cpp_min(lexmaxface[(i - 1) as usize], ub);
                    debug_assert!(sym_ge(
                        lexmaxface[(i - 1) as usize],
                        lexmaxface[i as usize],
                        eps,
                        kHighsInf
                    ));

                    // are we still equal?
                    if sym_gt(
                        lexmaxface[(i - 1) as usize],
                        lexmaxface[i as usize],
                        eps,
                        kHighsInf,
                    ) {
                        iseq = false;
                    } else if lexmaxepsrow[(colid - 1) as usize] == rowid {
                        debug_assert!(sym_eq(
                            lexmaxface[(i - 1) as usize],
                            lexmaxface[i as usize],
                            eps,
                            kHighsInf
                        ));
                        debug_assert!(
                            !is_integral
                                [orbidata.vars[origrowid as usize][origcolid as usize] as usize]
                        );
                        debug_assert!(!is_integral[var as usize]);
                        // left column (colid-1) has value x - epsilon, right column (colid) has value x, now
                        // must also become  x - epsilon in order to be larger or equal
                        // by axioms, we can squeeze infinitesimals between one other; epsilon > epsilon.
                        iseq = false;
                        debug_assert!(lexmaxepsrow[colid as usize] == -1);
                        lexmaxepsrow[colid as usize] = rowid;
                    }

                    // is there room left to decrease this variable further?
                    if is_integral[var as usize] {
                        // discrete type with unit steps: Remove one from the lexmax-value.
                        debug_assert!(
                            lexmaxface[i as usize] == (lexmaxface[i as usize] as i64) as f64
                        );
                        if sym_ge(lexmaxface[i as usize] - 1.0, lb, eps, kHighsInf) {
                            lastunfixed = rowid;
                        }
                    } else {
                        // continuous type: if we can subtract an infinitesimal value to the current lexmaxface[i] value,
                        // mark row as 'lastunfixed'
                        if sym_gt(lexmaxface[i as usize], lb, eps, kHighsInf) {
                            lastunfixed = rowid;
                        }
                    }
                }
            } else {
                // there had been a row before which breaks the equality-condition, choose maximally possible value
                lexmaxface[i as usize] = col_upper[var as usize];
            }

            rowid += 1;
            i += ncols;
        }

        colid += 1;
    }

    Some((lexmaxface, lexmaxepsrow))
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::propagateStaticOrbitope (lines 1884-1983), face comparison and fixing
pub(crate) fn static_orbitope_fix(
    ms: &mut HighsMipSolver,
    localdom: &mut HighsDomain,
    orbidata: &S3OrbitopeData,
    roworder: Option<&[HighsInt]>,
    nselrows: HighsInt,
    colorder: Option<&[HighsInt]>,
    is_integral: &[bool],
    eps: f64,
    feastol: f64,
) -> (bool, HighsInt) {
    if nselrows <= 0 {
        return (false, 0);
    }

    let ncols = orbidata.ncols;

    let faces = match static_orbitope_faces(
        orbidata,
        roworder,
        nselrows,
        colorder,
        &localdom.col_lower_,
        &localdom.col_upper_,
        is_integral,
        eps,
    ) {
        Some(f) => f,
        None => return (true, 0),
    };

    let lexminface = &faces.lexminface;
    let lexmaxface = &faces.lexmaxface;

    let mut nfixedvars: HighsInt = 0;

    for colid in 0..ncols {
        let mut i = colid;
        for rowid in 0..nselrows {
            debug_assert!(i == rowid * ncols + colid);

            let origrowid = get_array_entry_or_index(roworder, rowid);
            let origcolid = get_array_entry_or_index(colorder, colid);
            let var = orbidata.vars[origrowid as usize][origcolid as usize];
            let is_int = is_integral[var as usize];

            if sym_eq(
                lexminface[i as usize],
                lexmaxface[i as usize],
                eps,
                kHighsInf,
            ) {
                // tighten LB and UB to same value (i.e. fixing)
                let (infeasible, success) =
                    s3_tighten_lb(ms, localdom, var, lexminface[i as usize], is_int, feastol);
                if success {
                    nfixedvars += 1;
                }
                if infeasible {
                    return (true, nfixedvars);
                }

                let (infeasible, success) =
                    s3_tighten_ub(ms, localdom, var, lexminface[i as usize], is_int, feastol);
                if success {
                    nfixedvars += 1;
                }
                if infeasible {
                    return (true, nfixedvars);
                }
            } else {
                // lower bound, based on lexminface
                let (infeasible, success) =
                    s3_tighten_lb(ms, localdom, var, lexminface[i as usize], is_int, feastol);
                if success {
                    nfixedvars += 1;
                }
                if infeasible {
                    return (true, nfixedvars);
                }

                // upper bound, based on lexmaxface
                let (infeasible, success) =
                    s3_tighten_ub(ms, localdom, var, lexmaxface[i as usize], is_int, feastol);
                if success {
                    nfixedvars += 1;
                }
                if infeasible {
                    return (true, nfixedvars);
                }

                break;
            }

            i += ncols;
        }
    }

    (false, nfixedvars)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nine_col_solver() -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 9;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0; 9];
        lp.col_lower_ = vec![0.0; 9];
        lp.col_upper_ = vec![1.0; 9];
        lp.integrality_ = vec![HighsVarType::kInteger; 9];
        lp.a_matrix_.num_col_ = 9;
        lp.a_matrix_.num_row_ = 0;
        lp.a_matrix_.start_ = vec![0; 10];
        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        ms.mipdata_.cliquetable = crate::mip::Slot::new(HighsCliqueTable::new(9));
        let implications = HighsImplications::new(&ms);
        ms.mipdata_.implications = crate::mip::Slot::new(implications);
        ms
    }

    #[test]
    fn s3_face2_fixes_forced_column_ubs_3x3() {
        let mut ms = nine_col_solver();

        let mut localdom = HighsDomain::default();
        localdom.col_lower_ = vec![0.0; 9];
        localdom.col_upper_ = vec![1.0; 9];
        localdom.col_lower_pos_ = vec![-1; 9];
        localdom.col_upper_pos_ = vec![-1; 9];
        localdom.changedcolsflags_ = vec![0; 9];

        // Fix var 0 to 1 (lower bound)
        localdom.col_lower_[0] = 1.0;
        // Fix vars 1,2,4,7 to 0 (upper bound)
        localdom.col_upper_[1] = 0.0;
        localdom.col_upper_[2] = 0.0;
        localdom.col_upper_[4] = 0.0;
        localdom.col_upper_[7] = 0.0;

        let orbidata = S3OrbitopeData {
            vars: vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7, 8]],
            nrows: 3,
            ncols: 3,
            ..Default::default()
        };

        let is_integral = vec![true; 9];

        let (infeasible, nfixed) = static_orbitope_fix(
            &mut ms,
            &mut localdom,
            &orbidata,
            None,
            3,
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
    fn s3_face1_free_binary_2x2() {
        let orbidata = S3OrbitopeData {
            nrows: 2,
            ncols: 2,
            vars: vec![vec![0, 1], vec![2, 3]],
            ..Default::default()
        };
        let col_lower = vec![0.0f64; 4];
        let col_upper = vec![1.0f64; 4];
        let is_integral = vec![true; 4];
        let faces = static_orbitope_faces(
            &orbidata,
            None,
            2,
            None,
            &col_lower,
            &col_upper,
            &is_integral,
            1e-9,
        )
        .expect("expected feasible");
        assert_eq!(faces.lexminface, vec![0.0, 0.0, 0.0, 0.0]);
        assert_eq!(faces.lexmaxface, vec![1.0, 1.0, 1.0, 1.0]);
        assert_eq!(faces.lexminepsrow, vec![-1, -1]);
        assert_eq!(faces.lexmaxepsrow, vec![-1, -1]);
    }

    #[test]
    fn s3_face1_fixed_zero_entry_2x2() {
        let orbidata = S3OrbitopeData {
            nrows: 2,
            ncols: 2,
            vars: vec![vec![0, 1], vec![2, 3]],
            ..Default::default()
        };
        let col_lower = vec![0.0f64; 4];
        let mut col_upper = vec![1.0f64; 4];
        col_upper[0] = 0.0;
        let is_integral = vec![true; 4];
        let faces = static_orbitope_faces(
            &orbidata,
            None,
            2,
            None,
            &col_lower,
            &col_upper,
            &is_integral,
            1e-9,
        )
        .expect("expected feasible");
        assert_eq!(faces.lexminface, vec![0.0, 0.0, 0.0, 0.0]);
        assert_eq!(faces.lexmaxface, vec![0.0, 0.0, 1.0, 1.0]);
        assert_eq!(faces.lexminepsrow, vec![-1, -1]);
        assert_eq!(faces.lexmaxepsrow, vec![-1, -1]);
    }

    /// Exercises the repair path (lines 1611-1641 and 1766-1797) in both faces: row 1 forces
    /// a strict inequality that is repaired at row 0, then row 1 restarts unconstrained.
    #[test]
    fn s3_face1_repair_2x2() {
        let orbidata = S3OrbitopeData {
            nrows: 2,
            ncols: 2,
            vars: vec![vec![0, 1], vec![2, 3]],
            ..Default::default()
        };
        // x2 fixed to 0, x3 fixed to 1, x0 and x1 free binaries
        let col_lower = vec![0.0f64, 0.0, 0.0, 1.0];
        let col_upper = vec![1.0f64, 1.0, 0.0, 1.0];
        let is_integral = vec![true; 4];
        let faces = static_orbitope_faces(
            &orbidata,
            None,
            2,
            None,
            &col_lower,
            &col_upper,
            &is_integral,
            1e-9,
        )
        .expect("expected feasible");
        assert_eq!(faces.lexminface, vec![1.0, 0.0, 0.0, 1.0]);
        assert_eq!(faces.lexmaxface, vec![1.0, 0.0, 0.0, 1.0]);
        assert_eq!(faces.lexminepsrow, vec![-1, -1]);
        assert_eq!(faces.lexmaxepsrow, vec![-1, -1]);
    }

    #[test]
    fn s3_face1_infeasible_1x2() {
        let orbidata = S3OrbitopeData {
            nrows: 1,
            ncols: 2,
            vars: vec![vec![0, 1]],
            ..Default::default()
        };
        let col_lower = vec![0.0f64, 1.0];
        let col_upper = vec![0.0f64, 1.0];
        let is_integral = vec![true, true];
        let faces = static_orbitope_faces(
            &orbidata,
            None,
            1,
            None,
            &col_lower,
            &col_upper,
            &is_integral,
            1e-9,
        );
        assert!(faces.is_none());
    }
}
