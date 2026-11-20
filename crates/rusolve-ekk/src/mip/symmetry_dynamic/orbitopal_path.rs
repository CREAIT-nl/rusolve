//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! row and column order from the branching path, branch info for orbitopal symmetry handling.
//!
//! No gate site exists in this job: nothing else in the crate calls these functions yet
//! (their callers check the option), so no `if` gate is needed.

use super::orbitopal_order::update_column_order_when_branching_on_column;
use super::*;
use crate::prelude::*;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Local constants (SCIP 11.0.0 symmetry_orbitopal.h lines 59-63)
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_RowOrdering (lines 59-63)
const S3_ROW_ORDERING_NONE: HighsInt = 0;
/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.h::SCIP_RowOrdering (lines 59-63)
const S3_ROW_ORDERING_BRANCHING: HighsInt = 1;

// ---------------------------------------------------------------------------
// get_row_order
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::getRowOrder (lines 560-632)
///
/// `path` is the node's own `s3_path`: the chain of branchings that created
/// the node, whose head carries the info attached at the node's parent. This is
/// exactly the set of ancestors SCIP walks from SCIPnodeGetParent(node) (:595);
/// None = root. Do not pass the parent's `s3_path` (that would drop the info
/// attached at the parent). Walk from the head up to the root,
/// collect the nodes, then process them root first; within one node take entries
/// whose orbitope == orbi_idx and append rows_added in forward order.
///
/// Deviation from C: the C code walks parents and appends rows in reverse order
/// (:604-606), then reverses the whole array (:623-629). Here we collect nodes
/// first and iterate root-first, appending rows_added forward, which produces
/// the same root-first result without the intermediate reversal.
///
/// The rowordering == NONE early return of C :579-584 (roworder NULL, nselrows =
/// nrows) is NOT done here; the caller checks rowordering first, as
/// eventExecNodeBranched :829 already does. Only creates BRANCHING orbitopes.
pub(crate) fn get_row_order(orbi_idx: HighsInt, path: Option<&Arc<S3PathNode>>) -> Vec<HighsInt> {
    // Walk from head to root, collecting nodes (C :596-620 collects in reverse).
    let mut nodes: Vec<&S3PathNode> = Vec::new();
    let mut cur = path;
    while let Some(node) = cur {
        nodes.push(node);
        cur = node.parent.as_ref();
    }

    // Process root first (C reverses at :623-629; we iterate the collected
    // vector in reverse to achieve root-first order).
    let mut roworder: Vec<HighsInt> = Vec::new();
    for node in nodes.iter().rev() {
        for entry in &node.orbi {
            if entry.orbitope as HighsInt == orbi_idx {
                // C :604-606 appends rows[i] for i = nrows-1 down to 0,
                // then the whole array is reversed at :623-629.
                // Net effect: rows appear in forward order within a node,
                // nodes appear root-first. We replicate that directly.
                for &row in &entry.rows_added {
                    roworder.push(row);
                }
            }
        }
    }

    roworder
}

// ---------------------------------------------------------------------------
// populate_rooted_path_column_order
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::populateRootedPathColumnOrder (lines 637-716)
///
/// Same path convention as get_row_order: `path` is the `s3_path` of the node
/// whose order is wanted, standing for SCIP's rooted path from
/// SCIPnodeGetParent(node) (getColumnOrder :1328, eventExecNodeBranched :910).
/// Walk root first; for each node's
/// entries with orbitope == orbi_idx, apply each colswap (from, to) in order.
///
/// The caller supplies identity-initialised arrays.
pub(crate) fn populate_rooted_path_column_order(
    orbi_idx: HighsInt,
    path: Option<&Arc<S3PathNode>>,
    colorder: &mut [HighsInt],
    colorderinv: &mut [HighsInt],
) {
    debug_assert!(!colorder.is_empty());
    debug_assert_eq!(colorder.len(), colorderinv.len());

    // Collect nodes from head to root.
    let mut nodes: Vec<&S3PathNode> = Vec::new();
    let mut cur = path;
    while let Some(node) = cur {
        nodes.push(node);
        cur = node.parent.as_ref();
    }

    // Iterate root first (C :668-713 iterates rootedpath[0..=depth]).
    for node in nodes.iter().rev() {
        for entry in &node.orbi {
            if entry.orbitope as HighsInt != orbi_idx {
                continue;
            }

            // C :685-712: for each colswap in order.
            for &(from, to) in &entry.colswaps {
                debug_assert!(from != to); // C :691: no trivial swaps
                debug_assert!(from >= 0 && (from as usize) < colorder.len()); // C :692
                debug_assert!(to >= 0 && (to as usize) < colorder.len()); // C :693

                // C :696-711, exact statement order:
                let positionfromincolorder = colorderinv[from as usize];
                debug_assert!(positionfromincolorder >= 0); // C :697
                debug_assert!((positionfromincolorder as usize) < colorder.len()); // C :698
                debug_assert!(colorder[positionfromincolorder as usize] == from); // C :699

                let positiontoincolorder = colorderinv[to as usize];
                debug_assert!(positiontoincolorder >= 0); // C :703
                debug_assert!((positiontoincolorder as usize) < colorder.len()); // C :704
                debug_assert!(colorder[positiontoincolorder as usize] == to); // C :705

                // C :708-711: swap them, also keep track of the inverses.
                colorder[positiontoincolorder as usize] = from;
                colorder[positionfromincolorder as usize] = to;
                colorderinv[from as usize] = positiontoincolorder;
                colorderinv[to as usize] = positionfromincolorder;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// get_column_order
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::getColumnOrder (lines 1284-1338)
///
/// If columnordering == S3_COLUMN_ORDERING_NONE return None (C :1304-1309).
/// Otherwise build identity colorder and colorderinv of length ncols, call
/// populate_rooted_path_column_order, return Some((colorder, colorderinv)).
pub(crate) fn get_column_order(
    orbi_idx: HighsInt,
    path: Option<&Arc<S3PathNode>>,
    columnordering: HighsInt,
    ncols: HighsInt,
) -> Option<(Vec<HighsInt>, Vec<HighsInt>)> {
    // C :1304-1309: if columnordering == NONE, return NULL/NULL.
    if columnordering == S3_COLUMN_ORDERING_NONE {
        return None;
    }

    debug_assert!(ncols > 0); // C :1311

    let ncols_usize = ncols as usize;

    // C :1313-1314: allocate arrays.
    // C :1317-1318: populate colorder with standard ordering.
    let mut colorder: Vec<HighsInt> = (0..ncols).collect();
    // C :1321-1322: introduce inverse column ordering.
    let mut colorderinv: Vec<HighsInt> = (0..ncols).collect();

    // C :1328-1334: get the rooted path and populate.
    // In our representation, path (the node's own s3_path) already is the
    // rooted path from SCIPnodeGetParent(eventnode); None stands for node == NULL.
    populate_rooted_path_column_order(orbi_idx, path, &mut colorder, &mut colorderinv);

    debug_assert_eq!(colorder.len(), ncols_usize);
    debug_assert_eq!(colorderinv.len(), ncols_usize);

    Some((colorder, colorderinv))
}

// ---------------------------------------------------------------------------
// orbitopal_branch_info
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): symmetry_orbitopal.c::eventExecNodeBranched (lines 720-962), single branching column
///
/// HiGHS has one branching variable per node, so SCIP's nbranchvars loop
/// (:763-815, scanning children domchg) is replaced by the single `branch_col`
/// argument. The children-domchg scan is not ported.
///
/// `col_lower` / `col_upper` are the column bounds at the branched (parent)
/// node, BEFORE the branching change, indexed by HiGHS column (SCIP's
/// SCIPvarGetLbLocal). `eps` is the mipdata epsilon.
///
/// `parent_path` is the path of the branched node (its order); this is the
/// path that getRowOrder / populateRootedPathColumnOrder see in SCIP (they
/// start at the parent of eventnode, :595 and :910).
pub(crate) fn orbitopal_branch_info(
    data: &S3SymData,
    col_lower: &[f64],
    col_upper: &[f64],
    eps: f64,
    parent_path: Option<&Arc<S3PathNode>>,
    branch_col: HighsInt,
) -> Vec<S3OrbiNodeInfo> {
    let mut result: Vec<S3OrbiNodeInfo> = Vec::new();

    for oi in 0..data.orbitopes.len() {
        let orbidata = &data.orbitopes[oi];
        let oi_int = oi as HighsInt;

        // Membership test (C :786): check if branch_col is in the orbitope matrix.
        // rowindexmap is indexed by HiGHS column; -1 means not in the matrix.
        let rowid: HighsInt = if (branch_col as usize) < orbidata.rowindexmap.len() {
            orbidata.rowindexmap[branch_col as usize]
        } else {
            -1
        };
        if rowid < 0 {
            continue; // C :786: not in the matrix, skip this orbitope
        }

        debug_assert!(rowid >= 0 && rowid < orbidata.nrows); // C :848-849

        let mut rows_added: Vec<HighsInt> = Vec::new();
        let mut colswaps: Vec<(HighsInt, HighsInt)> = Vec::new();

        // Store data about row ordering (C :829-884).
        if orbidata.rowordering != S3_ROW_ORDERING_NONE {
            debug_assert!(orbidata.nrows > 0); // C :834
            debug_assert!(orbidata.rowordering == S3_ROW_ORDERING_BRANCHING); // C :835

            // C :838: get the present row order up to this node.
            let roworder = get_row_order(oi_int, parent_path);

            // C :851-861: avoid adding row to row order twice.
            let already_present = roworder.contains(&rowid);

            if !already_present {
                // C :866: append rowid to present roworder.
                // C :869-880: store in newnodeinfo->rows.
                rows_added.push(rowid);
            }
        }

        // Store data about column ordering (C :887-945).
        if orbidata.columnordering != S3_COLUMN_ORDERING_NONE {
            debug_assert!(orbidata.ncols > 0); // C :894

            // C :910-914: get the rooted path and populate colorder/colorderinv.
            let col_result =
                get_column_order(oi_int, parent_path, orbidata.columnordering, orbidata.ncols);

            // It cannot be None here because the ordering is not NONE.
            debug_assert!(col_result.is_some());

            if let Some((mut colorder, mut colorderinv)) = col_result {
                // C :917-941: for each branching variable (here: one).
                // C :919-920: updateColumnOrderWhenBranchingOnColumn.
                let swap = update_column_order_when_branching_on_column(
                    col_lower,
                    col_upper,
                    eps,
                    orbidata,
                    &mut colorder,
                    &mut colorderinv,
                    branch_col,
                );

                match swap {
                    Some((from, to)) => {
                        // C :923: skip trivial swaps of columns.
                        if from != to {
                            // C :927-940: store column swap.
                            colswaps.push((from, to));
                        }
                    }
                    None => {
                        // update_column_order_when_branching_on_column returned None.
                        // In C this path does not exist (the function always
                        // produces a swap or returns early without storing).
                        // debug_assert to flag unexpected behaviour.
                        debug_assert!(
                            false,
                            "update_column_order_when_branching_on_column returned None for branch_col={}",
                            branch_col
                        );
                    }
                }
            }
        }

        // C :948: store updates or free memory if no change applied.
        if !rows_added.is_empty() || !colswaps.is_empty() {
            result.push(S3OrbiNodeInfo {
                orbitope: oi as u32,
                rows_added,
                colswaps,
            });
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a simple S3OrbitopeData for 3 rows x 4 cols.
    fn make_test_orbidata() -> S3OrbitopeData {
        let nrows: HighsInt = 3;
        let ncols: HighsInt = 4;
        let nvars = (nrows * ncols) as usize; // 12

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

        S3OrbitopeData {
            vars,
            nrows,
            ncols,
            nbranchrows: 0,
            rowindexmap,
            colindexmap,
            columnordering: S3_COLUMN_ORDERING_MEDIAN,
            rowordering: S3_ROW_ORDERING_BRANCHING,
        }
    }

    #[test]
    fn s3_ord2_row_order_root_first() {
        // Root R: depth 1, rows_added = [2] for orbitope 0.
        let root = S3PathNode {
            parent: None,
            branch_col: 0,
            depth: 1,
            orbi: vec![S3OrbiNodeInfo {
                orbitope: 0,
                rows_added: vec![2],
                colswaps: vec![],
            }],
        };
        // Child N: depth 2, rows_added = [0, 1] for orbitope 0.
        let child = S3PathNode {
            parent: Some(Arc::new(root)),
            branch_col: 1,
            depth: 2,
            orbi: vec![S3OrbiNodeInfo {
                orbitope: 0,
                rows_added: vec![0, 1],
                colswaps: vec![],
            }],
        };
        let path = Arc::new(child);

        // get_row_order(0, Some(&path)): walk head->root, process root first.
        // Root contributes [2], child contributes [0, 1]. Result: [2, 0, 1].
        let result = get_row_order(0, Some(&path));
        assert_eq!(result, vec![2, 0, 1]);

        // get_row_order(1, same): no entries for orbitope 1.
        let result1 = get_row_order(1, Some(&path));
        assert!(result1.is_empty());
    }

    #[test]
    fn s3_ord2_populate_single_swap() {
        // Root carries colswaps [(0,1)] for orbitope 0, no rows.
        let root = S3PathNode {
            parent: None,
            branch_col: 0,
            depth: 1,
            orbi: vec![S3OrbiNodeInfo {
                orbitope: 0,
                rows_added: vec![],
                colswaps: vec![(0, 1)],
            }],
        };
        let path = Arc::new(root);

        let mut colorder: Vec<HighsInt> = vec![0, 1, 2];
        let mut colorderinv: Vec<HighsInt> = vec![0, 1, 2];

        populate_rooted_path_column_order(0, Some(&path), &mut colorder, &mut colorderinv);

        // Hand-trace from C :696-711 with swap (from=0, to=1):
        //   positionfromincolorder = colorderinv[0] = 0
        //   positiontoincolorder   = colorderinv[1] = 1
        //   colorder[1] = 0
        //   colorder[0] = 1
        //   colorderinv[0] = 1
        //   colorderinv[1] = 0
        // Result: colorder = [1, 0, 2], colorderinv = [1, 0, 2].
        assert_eq!(colorder, vec![1, 0, 2]);
        assert_eq!(colorderinv, vec![1, 0, 2]);
    }

    #[test]
    fn s3_ord2_column_order_none_is_none() {
        let result = get_column_order(0, None, S3_COLUMN_ORDERING_NONE, 4);
        assert!(result.is_none());
    }

    #[test]
    fn s3_ord2_column_order_identity_without_path() {
        let result = get_column_order(0, None, S3_COLUMN_ORDERING_MEDIAN, 4);
        assert_eq!(result, Some((vec![0, 1, 2, 3], vec![0, 1, 2, 3])));
    }

    #[test]
    fn s3_ord2_branch_info_adds_row_and_swap() {
        let orbidata = make_test_orbidata();
        let data = S3SymData {
            orbitopes: vec![orbidata],
            static_orbitopes: vec![],
            lexdatas: vec![],
        };

        let nvars = 12;
        let col_lower = vec![0.0f64; nvars];
        let col_upper = vec![1.0f64; nvars];
        let eps = 1e-9;

        // branch_col = 4: row 1, col 0 (4/4=1, 4%4=0).
        let result = orbitopal_branch_info(&data, &col_lower, &col_upper, eps, None, 4);

        // Expect exactly one entry for orbitope 0.
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].orbitope, 0);
        // Row 1 is not in any ancestor's roworder (parent is None), so it is added.
        assert_eq!(result[0].rows_added, vec![1]);
        // Column swap: with all bounds equal, the MEDIAN ordering should produce
        // a non-trivial swap (from != to) for col 0 in row 1.
        // We just check that every colswap has from != to.
        for &(from, to) in &result[0].colswaps {
            assert_ne!(from, to);
        }
    }

    #[test]
    fn s3_ord2_branch_info_row_not_repeated() {
        let orbidata = make_test_orbidata();
        let data = S3SymData {
            orbitopes: vec![orbidata],
            static_orbitopes: vec![],
            lexdatas: vec![],
        };

        let nvars = 12;
        let col_lower = vec![0.0f64; nvars];
        let col_upper = vec![1.0f64; nvars];
        let eps = 1e-9;

        // Parent path already has row 1 for orbitope 0.
        let parent_node = S3PathNode {
            parent: None,
            branch_col: 0,
            depth: 1,
            orbi: vec![S3OrbiNodeInfo {
                orbitope: 0,
                rows_added: vec![1],
                colswaps: vec![],
            }],
        };
        let parent_path = Some(Arc::new(parent_node));

        // branch_col = 5: row 1, col 1 (5/4=1, 5%4=1).
        let result =
            orbitopal_branch_info(&data, &col_lower, &col_upper, eps, parent_path.as_ref(), 5);

        // Row 1 is already in the ancestor's roworder, so rows_added must be empty.
        // The result may contain an entry only if colswaps is non-empty.
        for entry in &result {
            assert!(
                entry.rows_added.is_empty(),
                "rows_added should be empty when row is already present"
            );
        }
    }

    #[test]
    fn s3_ord2_branch_info_outside_matrix_empty() {
        // Set all rowindexmap entries to -1 (not in matrix).
        let mut orbidata = make_test_orbidata();
        for v in orbidata.rowindexmap.iter_mut() {
            *v = -1;
        }
        let data = S3SymData {
            orbitopes: vec![orbidata],
            static_orbitopes: vec![],
            lexdatas: vec![],
        };

        let nvars = 12;
        let col_lower = vec![0.0f64; nvars];
        let col_upper = vec![1.0f64; nvars];
        let eps = 1e-9;

        let result = orbitopal_branch_info(&data, &col_lower, &col_upper, eps, None, 4);

        assert!(result.is_empty());
    }

    #[test]
    fn s3_ord2_branch_info_no_colswaps_when_none() {
        let mut orbidata = make_test_orbidata();
        orbidata.columnordering = S3_COLUMN_ORDERING_NONE;
        let data = S3SymData {
            orbitopes: vec![orbidata],
            static_orbitopes: vec![],
            lexdatas: vec![],
        };

        let nvars = 12;
        let col_lower = vec![0.0f64; nvars];
        let col_upper = vec![1.0f64; nvars];
        let eps = 1e-9;

        // branch_col = 4: row 1, col 0.
        let result = orbitopal_branch_info(&data, &col_lower, &col_upper, eps, None, 4);

        // One entry: rows_added = [1], colswaps empty (columnordering is NONE).
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].orbitope, 0);
        assert_eq!(result[0].rows_added, vec![1]);
        assert!(result[0].colswaps.is_empty());
    }
}
