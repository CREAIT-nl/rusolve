//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! lexdataCreate, variable order, and the path-based replacement of
//! shadowtreeFillNodeDepthBranchIndices for the dynamic symmetry propagation
//! behind `mip_symmetry_dynamic`.

use crate::mip::symmetry_dynamic::{S3LexData, S3PathNode};
use crate::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// NodeDepthBranchIndex
// ---------------------------------------------------------------------------

/// (depth, index) information for a variable on the rooted branching path.
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::NodeDepthBranchIndex (lines 120-125)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct NodeDepthBranchIndex {
    pub(crate) nodedepth: HighsInt,
    pub(crate) index: HighsInt,
}

// ---------------------------------------------------------------------------
// lexdata_create
// ---------------------------------------------------------------------------

/// Port of lexdataCreate: build an S3LexData from a permutation, removing
/// fixed points. Unsigned permutations only (Q7: no signed branch).
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::lexdataCreate (lines 222-405)
pub(crate) fn lexdata_create(
    symvarmap: &mut HashMap<HighsInt, HighsInt>,
    nsymvars: &mut HighsInt,
    vars: &[HighsInt],
    perm: &[HighsInt],
    usedynamicorder: bool,
) -> Option<S3LexData> {
    let nvars = vars.len();
    debug_assert!(perm.len() == nvars);

    // Remove fixed points: indexcorrection (C lines 264-272).
    let mut indexcorrection = vec![0i32; nvars];
    let mut naffectedvariables: HighsInt = 0;
    for i in 0..nvars {
        if perm[i] == i as HighsInt {
            indexcorrection[i] = -1;
        } else {
            indexcorrection[i] = naffectedvariables;
            naffectedvariables += 1;
        }
    }

    // Do nothing if reductions would be trivial (C lines 275-283).
    if naffectedvariables <= 0 {
        debug_assert!(naffectedvariables == 0);
        return None;
    }

    // C lines 286-291: SCIPactivateShadowTree is NOT ported; the path walk
    // (lexred_branch_order_from_path) replaces the shadow tree.
    // C lines 363-369: SCIPmarkDoNotMultaggrVar and SCIPcaptureVar are dropped
    // (no SCIP memory management).

    let nvars_affected = naffectedvariables as usize;

    let mut lexdata = S3LexData {
        isdynamic: usedynamicorder,
        vars: vec![0i32; nvars_affected],
        perm: vec![0i32; nvars_affected],
        invperm: vec![0i32; nvars_affected],
        varmap: HashMap::new(),
        symvarmap: Vec::new(),
        nsymvars: 0,
        hasdynamicperm: usedynamicorder,
    };

    // Determine the vars, perm (C lines 308-341). Unsigned only: the
    // issignedperm branches are omitted.
    for j in 0..nvars {
        let i = indexcorrection[j];
        if i < 0 {
            continue;
        }
        let iu = i as usize;

        // j is the original index, i is the relabeled index.
        lexdata.vars[iu] = vars[j];
        lexdata.perm[iu] = indexcorrection[perm[j] as usize];

        debug_assert!(perm[j] != j as HighsInt);
        debug_assert!(lexdata.perm[iu] != i);
        debug_assert!(lexdata.perm[iu] >= 0);
        debug_assert!(lexdata.perm[iu] < nvars_affected as HighsInt);
    }

    // Determine invperm (C lines 344-360). Unsigned only.
    for i in 0..nvars_affected {
        lexdata.invperm[lexdata.perm[i] as usize] = i as HighsInt;
    }

    // Create hashmap for all variables (C lines 371-402).
    if usedynamicorder {
        for i in 0..nvars_affected {
            let var = lexdata.vars[i];

            // The hashmap in lexdata maps to the index in the vars array.
            lexdata.varmap.insert(var, i as HighsInt);

            // Var already added to symvarmap?
            if symvarmap.contains_key(&var) {
                continue;
            }

            // The hashmap in symvarmap maps to a unique index.
            symvarmap.insert(var, *nsymvars);
            *nsymvars += 1;
        }
    }
    // else: varmap stays empty (C line 402: (*lexdata)->varmap = NULL).

    Some(lexdata)
}

// ---------------------------------------------------------------------------
// sort_by_node_depth_branch_indices
// ---------------------------------------------------------------------------

/// Comparator: by nodedepth ascending, then index ascending.
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::sortbynodedepthbranchindices (lines 408-461)
pub(crate) fn sort_by_node_depth_branch_indices(
    index1: &NodeDepthBranchIndex,
    index2: &NodeDepthBranchIndex,
) -> Ordering {
    if index1.nodedepth < index2.nodedepth {
        Ordering::Less
    } else if index1.nodedepth > index2.nodedepth {
        Ordering::Greater
    } else {
        debug_assert!(index1.index != index2.index);
        if index1.index < index2.index {
            Ordering::Less
        } else if index1.index > index2.index {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

// ---------------------------------------------------------------------------
// var_order_get_index
// ---------------------------------------------------------------------------

/// Return the index of a variable at a specific position of a variable order.
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::varOrderGetIndex (lines 464-476)
pub(crate) fn var_order_get_index(varorder: Option<&[HighsInt]>, pos: HighsInt) -> HighsInt {
    debug_assert!(pos >= 0);
    match varorder {
        None => pos,
        Some(order) => order[pos as usize],
    }
}

// ---------------------------------------------------------------------------
// lexred_branch_order_from_path
// ---------------------------------------------------------------------------

/// Walk the branching path from the focus node toward the root, recording the
/// shallowest (depth, index) for each branched column. Replaces
/// shadowtreeFillNodeDepthBranchIndices.
///
/// Unlike C (:1712-1716), columns outside symvarmap are not skipped, so
/// branchvars may be longer than C's. This only changes which branch of
/// get_var_order is taken; both branches select the same set (columns of the
/// permutation with nodedepth > 0), and the sort keys are unique, so the
/// resulting order is identical.
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::shadowtreeFillNodeDepthBranchIndices (lines 1636-1751)
pub(crate) fn lexred_branch_order_from_path(
    path: Option<&Arc<S3PathNode>>,
) -> (HashMap<HighsInt, NodeDepthBranchIndex>, Vec<HighsInt>) {
    let mut map: HashMap<HighsInt, NodeDepthBranchIndex> = HashMap::new();
    let mut branchvars: Vec<HighsInt> = Vec::new();

    let mut node = path;
    while let Some(n) = node {
        let col = n.branch_col;
        let d = n.depth;
        debug_assert!(d > 0);

        // First seen (leaf to root): add to branchvars.
        if !map.contains_key(&col) {
            branchvars.push(col);
        }

        // Overwrite: the shallowest depth wins (C: nodedepth == 0 ||
        // nodedepth > shadowdepth). index is always 0 because HiGHS branches
        // on one variable per node.
        map.insert(
            col,
            NodeDepthBranchIndex {
                nodedepth: d,
                index: 0,
            },
        );

        node = n.parent.as_ref();
    }

    (map, branchvars)
}

// ---------------------------------------------------------------------------
// get_var_order
// ---------------------------------------------------------------------------

/// Build the variable order for a LexRed permutation from the branching
/// information on the path.
///
/// SCIP 11.0.0 (32dd425): symmetry_lexred.c::getVarOrder (lines 481-585)
pub(crate) fn get_var_order(
    lexdata: &S3LexData,
    nodedepthbranchindices: &HashMap<HighsInt, NodeDepthBranchIndex>,
    branchvars: &[HighsInt],
    varorder: &mut Vec<HighsInt>,
) {
    // C line 516: *nselvars = 0.
    varorder.clear();

    let nvars = lexdata.vars.len();
    let nbranchvars = branchvars.len();

    if nvars < nbranchvars {
        // C lines 518-542: for permutations with small support, check all
        // support entries.
        for i in 0..nvars {
            let var = lexdata.vars[i];
            let nodedepth = nodedepthbranchindices.get(&var).map_or(0, |e| e.nodedepth);

            // Skip variables that have not been used for branching.
            if nodedepth <= 0 {
                continue;
            }

            varorder.push(i as HighsInt);
        }
    } else {
        // C lines 544-572: for permutations where the support is larger than
        // the number of branched vars, check the branched vars.
        for &var in branchvars {
            // C line 561: get the variable index in the lexdata->vars array.
            // C line 565: skip if not in the permutation (INT_MAX).
            if let Some(&varindex) = lexdata.varmap.get(&var) {
                varorder.push(varindex);
            }
        }
    }

    // C lines 574-582: sort the first n elements of varorder by depth, then
    // by index. Keys are unique by construction (each depth has at most one
    // branching variable on a HiGHS path), so the result is the same as
    // SCIPsortInd for any correct sort. Use sort_by (stable) rather than
    // sort_unstable.
    if varorder.len() > 1 {
        varorder.sort_by(|&a, &b| {
            let nda = nodedepthbranchindices
                .get(&lexdata.vars[a as usize])
                .expect("selected variable has a (depth, index) entry");
            let ndb = nodedepthbranchindices
                .get(&lexdata.vars[b as usize])
                .expect("selected variable has a (depth, index) entry");
            sort_by_node_depth_branch_indices(nda, ndb)
        });
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s3_lxd_fixed_point_removal() {
        // vars = [10, 11, 12, 13], perm = [0, 2, 1, 3]: positions 0 and 3 are
        // fixed points; positions 1 and 2 swap.
        let vars = vec![10, 11, 12, 13];
        let perm = vec![0, 2, 1, 3];

        // usedynamicorder = false: varmap empty, symvarmap unchanged.
        let mut symvarmap = HashMap::new();
        let mut nsymvars: HighsInt = 0;
        let result = lexdata_create(&mut symvarmap, &mut nsymvars, &vars, &perm, false);
        let lexdata = result.expect("permutation has affected positions");
        assert_eq!(lexdata.vars, vec![11, 12]);
        // Relabelled perm: position 1 (col 11) maps to position 2 (col 12),
        // i.e. relabelled index 0 maps to 1 and 1 maps to 0.
        assert_eq!(lexdata.perm, vec![1, 0]);
        assert_eq!(lexdata.invperm, vec![1, 0]);
        assert!(lexdata.varmap.is_empty());
        assert_eq!(nsymvars, 0);

        // All-fixed case: returns None.
        let vars_fixed = vec![1, 2];
        let perm_fixed = vec![0, 1];
        let mut symvarmap2 = HashMap::new();
        let mut nsymvars2: HighsInt = 0;
        assert!(
            lexdata_create(
                &mut symvarmap2,
                &mut nsymvars2,
                &vars_fixed,
                &perm_fixed,
                false
            )
            .is_none()
        );

        // usedynamicorder = true: varmap and symvarmap populated.
        let mut symvarmap3 = HashMap::new();
        let mut nsymvars3: HighsInt = 0;
        let result3 = lexdata_create(&mut symvarmap3, &mut nsymvars3, &vars, &perm, true)
            .expect("permutation has affected positions");
        assert_eq!(result3.varmap.get(&11), Some(&0));
        assert_eq!(result3.varmap.get(&12), Some(&1));
        assert_eq!(symvarmap3.get(&11), Some(&0));
        assert_eq!(symvarmap3.get(&12), Some(&1));
        assert_eq!(nsymvars3, 2);
    }

    #[test]
    fn s3_lxd_branch_order_support_smaller() {
        // Path: root depth 1 branch_col 12, then depth 2 branch_col 10, then
        // depth 3 branch_col 11. Head = depth-3 node.
        let root = Arc::new(S3PathNode {
            parent: None,
            branch_col: 12,
            depth: 1,
            orbi: Vec::new(),
        });
        let mid = Arc::new(S3PathNode {
            parent: Some(root.clone()),
            branch_col: 10,
            depth: 2,
            orbi: Vec::new(),
        });
        let head = Arc::new(S3PathNode {
            parent: Some(mid.clone()),
            branch_col: 11,
            depth: 3,
            orbi: Vec::new(),
        });

        let (map, branchvars) = lexred_branch_order_from_path(Some(&head));

        // Depths: 12->1, 10->2, 11->3, all index 0.
        assert_eq!(
            map[&12],
            NodeDepthBranchIndex {
                nodedepth: 1,
                index: 0
            }
        );
        assert_eq!(
            map[&10],
            NodeDepthBranchIndex {
                nodedepth: 2,
                index: 0
            }
        );
        assert_eq!(
            map[&11],
            NodeDepthBranchIndex {
                nodedepth: 3,
                index: 0
            }
        );

        // branchvars order: first seen leaf to root = [11, 10, 12].
        assert_eq!(branchvars, vec![11, 10, 12]);

        // LexData: vars [10, 11], perm [1, 0], usedynamicorder true.
        let mut symvarmap = HashMap::new();
        let mut nsymvars: HighsInt = 0;
        let lexdata = lexdata_create(&mut symvarmap, &mut nsymvars, &[10, 11], &[1, 0], true)
            .expect("permutation has affected positions");

        // nvars = 2 < nbranchvars = 3, so branch A.
        // i=0: var=10, nodedepth=2 > 0, push 0.
        // i=1: var=11, nodedepth=3 > 0, push 1.
        // Sort by depth ascending: 10 (depth 2) before 11 (depth 3).
        // Result: [0, 1].
        let mut varorder = Vec::new();
        get_var_order(&lexdata, &map, &branchvars, &mut varorder);
        assert_eq!(varorder, vec![0, 1]);
    }

    #[test]
    fn s3_lxd_branch_order_support_larger() {
        // Support of 4 columns [20, 21, 22, 23] with perm [1, 0, 3, 2].
        // Path branches on 22 (depth 1) and 21 (depth 2). nbranchvars = 2 <
        // nvars = 4, so branch B.
        let root = Arc::new(S3PathNode {
            parent: None,
            branch_col: 22,
            depth: 1,
            orbi: Vec::new(),
        });
        let head = Arc::new(S3PathNode {
            parent: Some(root.clone()),
            branch_col: 21,
            depth: 2,
            orbi: Vec::new(),
        });

        let (map, branchvars) = lexred_branch_order_from_path(Some(&head));

        // LexData: vars [20, 21, 22, 23], perm [1, 0, 3, 2], usedynamicorder
        // true. varmap: 20->0, 21->1, 22->2, 23->3.
        let mut symvarmap = HashMap::new();
        let mut nsymvars: HighsInt = 0;
        let lexdata = lexdata_create(
            &mut symvarmap,
            &mut nsymvars,
            &[20, 21, 22, 23],
            &[1, 0, 3, 2],
            true,
        )
        .expect("permutation has affected positions");

        // Branch B: iterate branchvars = [21, 22].
        // var=21: varmap gives 1, push 1.
        // var=22: varmap gives 2, push 2.
        // varorder before sort: [1, 2].
        // Sort by depth: 21 has depth 2, 22 has depth 1. Ascending: 22 (depth
        // 1) before 21 (depth 2). Result: [2, 1].
        let mut varorder = Vec::new();
        get_var_order(&lexdata, &map, &branchvars, &mut varorder);
        assert_eq!(varorder, vec![2, 1]);
    }

    #[test]
    fn s3_lxd_path_none_is_empty() {
        let (map, branchvars) = lexred_branch_order_from_path(None);
        assert!(map.is_empty());
        assert!(branchvars.is_empty());
    }
}
