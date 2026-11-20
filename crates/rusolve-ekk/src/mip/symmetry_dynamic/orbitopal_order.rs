//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! column equivalence and column-order update for orbitopal symmetry handling.
//!
//! No gate site exists in this job: nothing else in the crate calls these functions yet,
//! so no `if` gate is needed.

use super::*;
use crate::prelude::*;

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::testColumnsAreSymmetricallyEquivalent (lines 240-275)
pub(crate) fn test_columns_are_symmetrically_equivalent(
    col_lower: &[f64],
    col_upper: &[f64],
    orbidata: &S3OrbitopeData,
    col1: HighsInt,
    col2: HighsInt,
    eps: f64,
) -> bool {
    debug_assert!(col1 >= 0);
    debug_assert!(col1 < orbidata.ncols);
    debug_assert!(col2 >= 0);
    debug_assert!(col2 < orbidata.ncols);

    // @todo test only for the selected rows (see function description)
    for i in 0..orbidata.nrows {
        let var1 = orbidata.vars[i as usize][col1 as usize];
        let var2 = orbidata.vars[i as usize][col2 as usize];

        // if variable bounds differ: columns c and origcolid are not the same
        if !sym_eq(
            col_lower[var1 as usize],
            col_lower[var2 as usize],
            eps,
            kHighsInf,
        ) || !sym_eq(
            col_upper[var1 as usize],
            col_upper[var2 as usize],
            eps,
            kHighsInf,
        ) {
            return false;
        }
    }

    // loop terminated, so columns are equal
    true
}

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::updateColumnOrderWhenBranchingOnColumn (lines 290-509)
pub(crate) fn update_column_order_when_branching_on_column(
    col_lower: &[f64],
    col_upper: &[f64],
    eps: f64,
    orbidata: &S3OrbitopeData,
    colorder: &mut [HighsInt],
    colorderinv: &mut [HighsInt],
    var: HighsInt,
) -> Option<(HighsInt, HighsInt)> {
    debug_assert!(orbidata.ncols > 0);
    debug_assert!(orbidata.nrows > 0);

    let ncols = orbidata.ncols;
    debug_assert!(ncols > 0);

    // do not apply a column swap if no column permutations are applied
    if orbidata.columnordering == S3_COLUMN_ORDERING_NONE {
        return Some((0, 0));
    }

    // only variables from the orbitope matrix are of interest
    let origcolid = match orbidata.colindexmap.get(var as usize) {
        Some(&v) if v >= 0 => v,
        _ => return Some((0, 0)),
    };
    debug_assert!(origcolid >= 0);
    debug_assert!(origcolid < ncols);

    // policy: swap with identical column that is closest to the center in relabeled order
    // @todo other policies: If the variable is in a ppc-row, then select the minimal/second to minimal to branch on
    let mut swaporigcolid = origcolid;

    match orbidata.columnordering {
        S3_COLUMN_ORDERING_CENTRE | S3_COLUMN_ORDERING_FIRST | S3_COLUMN_ORDERING_LAST => {
            let middlecolumn = if orbidata.columnordering == S3_COLUMN_ORDERING_CENTRE {
                ncols / 2
            } else {
                0
            };

            // for each column, test column ordering condition, then if the column is identical to the var-column
            for col_c in 0..ncols {
                // origcolid is not interesting
                if col_c == origcolid {
                    continue;
                }

                // test if c is a better choice than swaporigcolid,
                // otherwise continue to next iteration through CONDITIONFAIL
                match orbidata.columnordering {
                    S3_COLUMN_ORDERING_FIRST => {
                        // only swap with c if c is earlier in column order than swaporigcolid
                        if colorderinv[col_c as usize] >= colorderinv[swaporigcolid as usize] {
                            continue;
                        }
                    }
                    S3_COLUMN_ORDERING_LAST => {
                        // only swap with c if c is later in column order than swaporigcolid
                        if colorderinv[col_c as usize] <= colorderinv[swaporigcolid as usize] {
                            continue;
                        }
                    }
                    S3_COLUMN_ORDERING_CENTRE => {
                        // if the column is not more central than swaporigcolid, ignore
                        if (colorderinv[col_c as usize] - middlecolumn).abs()
                            >= (colorderinv[swaporigcolid as usize] - middlecolumn).abs()
                        {
                            continue;
                        }
                    }
                    _ => return None,
                }

                // test: are c and origcolid the same columns w.r.t. the variable domain restrictions?
                if !test_columns_are_symmetrically_equivalent(
                    col_lower, col_upper, orbidata, col_c, origcolid, eps,
                ) {
                    continue;
                }

                // the variable domain reductions in c and origcolid are the same
                swaporigcolid = col_c;
            }
        }
        S3_COLUMN_ORDERING_MEDIAN => {
            // collect columns identical to the var-column, then select column satisfying ordering condition
            let mut origequalcolids: Vec<HighsInt> = Vec::new();

            // collect equal columns
            for col_c in 0..ncols {
                // column origcolid is always equal to itself
                if col_c == origcolid {
                    origequalcolids.push(col_c);
                    debug_assert!(origequalcolids.len() <= ncols as usize);
                    continue;
                }

                // test: are c and origcolid the same columns w.r.t. the variable domain restrictions?
                if !test_columns_are_symmetrically_equivalent(
                    col_lower, col_upper, orbidata, col_c, origcolid, eps,
                ) {
                    continue;
                }

                // the variable domain reductions in c and origcolid are the same
                origequalcolids.push(col_c);
                debug_assert!(origequalcolids.len() <= ncols as usize);
            }

            // we should have found origcolid, at least
            debug_assert!(!origequalcolids.is_empty());

            // from origequalcolids, select the column satisfying the column ordering policy

            // get median column; since colorder maps origcolids to our ordering,
            // we need to use colorderinv as the argument.
            // @todo computing the median is O(n) by repeated median-of-medians (CLRS, Ch9), but let's just sort things
            pdqsort_by(&mut origequalcolids, |a, b| {
                colorderinv[*a as usize] < colorderinv[*b as usize]
            });
            // get the median, that is swaporigcolid
            swaporigcolid = origequalcolids[origequalcolids.len() / 2];
        }
        _ => {
            // unknown column ordering variant
            return None;
        }
    }

    let from = swaporigcolid;
    let to = origcolid;

    // if we do not replace origcolid
    if swaporigcolid == origcolid {
        return Some((from, to));
    }

    // NDEBUG-only blocks at lines 309-314, 327-330, 474-483 are dropped.

    // now swap the permuted column indices of swaporigcolid and origcolid

    // at which column is origcolid?
    let positionorigcolidincolorder = colorderinv[origcolid as usize];
    debug_assert!(positionorigcolidincolorder >= 0);
    debug_assert!(positionorigcolidincolorder < ncols);
    debug_assert!(colorder[positionorigcolidincolorder as usize] == origcolid);

    // at which column is swaporigcolid?
    let positionswaporigcolidincolorder = colorderinv[swaporigcolid as usize];
    debug_assert!(positionswaporigcolidincolorder >= 0);
    debug_assert!(positionswaporigcolidincolorder < ncols);
    debug_assert!(colorder[positionswaporigcolidincolorder as usize] == swaporigcolid);

    // swap them, also keep track of the inverses
    colorder[positionswaporigcolidincolorder as usize] = origcolid;
    colorder[positionorigcolidincolorder as usize] = swaporigcolid;
    colorderinv[origcolid as usize] = positionswaporigcolidincolorder;
    colorderinv[swaporigcolid as usize] = positionorigcolidincolorder;

    Some((from, to))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_orbidata(columnordering: HighsInt) -> S3OrbitopeData {
        let vars: Vec<Vec<HighsInt>> = (0..3)
            .map(|i| (0..4).map(|j| (i * 4 + j) as HighsInt).collect())
            .collect();
        let colindexmap: Vec<HighsInt> = (0..12).map(|v| v % 4).collect();
        S3OrbitopeData {
            vars,
            nrows: 3,
            ncols: 4,
            nbranchrows: 0,
            rowindexmap: vec![],
            colindexmap,
            columnordering,
            rowordering: 0,
        }
    }

    #[test]
    fn s3_ord1_test_columns_equivalent() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_MEDIAN);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];

        assert!(test_columns_are_symmetrically_equivalent(
            &col_lower, &col_upper, &orbidata, 0, 3, 1e-9
        ));

        let mut col_lower2 = vec![0.0; 12];
        col_lower2[orbidata.vars[1][1] as usize] = 1.0;
        assert!(!test_columns_are_symmetrically_equivalent(
            &col_lower2,
            &col_upper,
            &orbidata,
            0,
            1,
            1e-9
        ));
    }

    #[test]
    fn s3_ord1_median_all_equal() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_MEDIAN);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            0,
        );
        assert_eq!(result, Some((2, 0)));
        assert_eq!(colorder, vec![2, 1, 0, 3]);
        assert_eq!(colorderinv, vec![2, 1, 0, 3]);
    }

    #[test]
    fn s3_ord1_median_excludes_unequal() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_MEDIAN);
        let mut col_lower = vec![0.0; 12];
        col_lower[orbidata.vars[0][1] as usize] = 1.0;
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            0,
        );
        assert_eq!(result, Some((2, 0)));
        assert_eq!(colorder, vec![2, 1, 0, 3]);
        assert_eq!(colorderinv, vec![2, 1, 0, 3]);
    }

    #[test]
    fn s3_ord1_first_picks_leftmost() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_FIRST);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            2,
        );
        assert_eq!(result, Some((0, 2)));
        assert_eq!(colorder, vec![2, 1, 0, 3]);
        assert_eq!(colorderinv, vec![2, 1, 0, 3]);
    }

    #[test]
    fn s3_ord1_last_picks_rightmost() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_LAST);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            0,
        );
        assert_eq!(result, Some((3, 0)));
        assert_eq!(colorder, vec![3, 1, 2, 0]);
        assert_eq!(colorderinv, vec![3, 1, 2, 0]);
    }

    #[test]
    fn s3_ord1_none_is_noop() {
        let orbidata = make_orbidata(S3_COLUMN_ORDERING_NONE);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            0,
        );
        assert_eq!(result, Some((0, 0)));
        assert_eq!(colorder, vec![0, 1, 2, 3]);
        assert_eq!(colorderinv, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s3_ord1_not_in_orbitope_is_noop() {
        let mut orbidata = make_orbidata(S3_COLUMN_ORDERING_MEDIAN);
        orbidata.colindexmap[5] = -1;
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            5,
        );
        assert_eq!(result, Some((0, 0)));
        assert_eq!(colorder, vec![0, 1, 2, 3]);
        assert_eq!(colorderinv, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s3_ord1_unknown_ordering_is_error() {
        let orbidata = make_orbidata(99);
        let col_lower = vec![0.0; 12];
        let col_upper = vec![1.0; 12];
        let mut colorder = vec![0, 1, 2, 3];
        let mut colorderinv = vec![0, 1, 2, 3];

        let result = update_column_order_when_branching_on_column(
            &col_lower,
            &col_upper,
            1e-9,
            &orbidata,
            &mut colorder,
            &mut colorderinv,
            0,
        );
        assert_eq!(result, None);
    }
}
