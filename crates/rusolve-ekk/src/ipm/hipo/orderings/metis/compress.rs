//! HiGHS source: extern/metis/libmetis/compress.c, extern/metis/libmetis/contig.c, extern/metis/libmetis/separator.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Graph compression (identical-adjacency vertices, keyed by `ikvsorti` on
//! adjacency checksums), dense-column pruning, separator-induced components
//! and the edge-to-vertex separator conversion.

use super::gklib::{icopy, ikv_t, ikvsorti, iset};
use super::graph::{CreateGraph, FreeRData, SetupGraph_label, SetupGraph_tvwgt};
use super::refine::{Allocate2WayNodePartitionMemory, Compute2WayNodePartitionParams};
use super::sfm::{FM_2WayNodeRefine1Sided, FM_2WayNodeRefine2Sided};
use super::{COMPRESSION_FRACTION, ctrl_t, graph_t, idx_t, real_t};

/// HiGHS: compress.c::CompressGraph (lines 25-141)
/// Returns `None` when compression would not reduce the graph below
/// `COMPRESSION_FRACTION` (the C `NULL`).
pub fn CompressGraph(
    ctrl: &mut ctrl_t,
    nvtxs: idx_t,
    xadj: &[idx_t],
    adjncy: &[idx_t],
    vwgt: Option<&[idx_t]>,
    cptr: &mut [idx_t],
    cind: &mut [idx_t],
) -> Option<graph_t> {
    // IFSET(ctrl->dbglvl, METIS_DBG_INFO, ...) (compress.c 82-83) dropped

    let mut mark = vec![-1 as idx_t; nvtxs as usize];
    let mut map = vec![-1 as idx_t; nvtxs as usize];
    let mut keys: Vec<ikv_t> = vec![ikv_t { key: 0, val: 0 }; nvtxs as usize];

    // Compute a key for each adjacency list
    for i in 0..nvtxs {
        let mut k: idx_t = 0;
        let mut j = xadj[i as usize];
        while j < xadj[(i + 1) as usize] {
            k += adjncy[j as usize];
            j += 1;
        }
        keys[i as usize].key = k + i; /* Add the diagonal entry as well */
        keys[i as usize].val = i;
    }

    ikvsorti(nvtxs as usize, &mut keys);

    let mut l: idx_t = 0;
    cptr[0] = 0;
    let mut cnvtxs: idx_t = 0;
    let mut i: idx_t = 0;
    while i < nvtxs {
        let ii = keys[i as usize].val;
        if map[ii as usize] == -1 {
            mark[ii as usize] = i; /* Add the diagonal entry */
            let mut j = xadj[ii as usize];
            while j < xadj[(ii + 1) as usize] {
                mark[adjncy[j as usize] as usize] = i;
                j += 1;
            }

            map[ii as usize] = cnvtxs;
            cind[l as usize] = ii;
            l += 1;

            let mut j = i + 1;
            while j < nvtxs {
                let iii = keys[j as usize].val;

                if keys[i as usize].key != keys[j as usize].key
                    || xadj[(ii + 1) as usize] - xadj[ii as usize]
                        != xadj[(iii + 1) as usize] - xadj[iii as usize]
                {
                    break; /* Break if keys or degrees are different */
                }

                if map[iii as usize] == -1 {
                    /* Do a comparison if iii has not been mapped */
                    let mut jj = xadj[iii as usize];
                    while jj < xadj[(iii + 1) as usize] {
                        if mark[adjncy[jj as usize] as usize] != i {
                            break;
                        }
                        jj += 1;
                    }

                    if jj == xadj[(iii + 1) as usize] {
                        /* Identical adjacency structure */
                        map[iii as usize] = cnvtxs;
                        cind[l as usize] = iii;
                        l += 1;
                    }
                }
                j += 1;
            }

            cnvtxs += 1;
            cptr[cnvtxs as usize] = l;
        }
        i += 1;
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_INFO, ...) (compress.c 82-83) dropped

    let mut graph: Option<graph_t> = None;

    // C: cnvtxs < COMPRESSION_FRACTION*nvtxs
    // COMPRESSION_FRACTION is f64, nvtxs is idx_t (i32).
    // C promotes nvtxs to f64, multiplies in f64, compares in f64.
    if (cnvtxs as f64) < COMPRESSION_FRACTION * (nvtxs as f64) {
        /* Sufficient compression is possible, so go ahead and create the
        compressed graph */

        let mut g = CreateGraph();

        let mut cnedges: idx_t = 0;
        for i in 0..cnvtxs {
            let ii = cind[cptr[i as usize] as usize];
            cnedges += xadj[(ii + 1) as usize] - xadj[ii as usize];
        }

        /* Allocate memory for the compressed graph */
        g.xadj = vec![0 as idx_t; (cnvtxs + 1) as usize];
        g.vwgt = vec![0 as idx_t; cnvtxs as usize];
        g.adjncy = vec![0 as idx_t; cnedges as usize];
        g.adjwgt = vec![1 as idx_t; cnedges as usize];

        /* Now go and compress the graph */
        iset(nvtxs as usize, -1, &mut mark);
        let mut l: idx_t = 0;
        g.xadj[0] = 0;
        for i in 0..cnvtxs {
            mark[i as usize] = i; /* Remove any diagonal entries in the compressed graph */
            let mut j = cptr[i as usize];
            while j < cptr[(i + 1) as usize] {
                let ii = cind[j as usize];

                /* accumulate the vertex weights of the constituent vertices */
                // C: cvwgt[i] += (vwgt == NULL ? 1 : vwgt[ii]);
                g.vwgt[i as usize] += match vwgt {
                    None => 1,
                    Some(v) => v[ii as usize],
                };

                /* generate the combined adjacency list */
                let mut jj = xadj[ii as usize];
                while jj < xadj[(ii + 1) as usize] {
                    let k = map[adjncy[jj as usize] as usize];
                    if mark[k as usize] != i {
                        mark[k as usize] = i;
                        g.adjncy[l as usize] = k;
                        l += 1;
                    }
                    jj += 1;
                }
                j += 1;
            }
            g.xadj[(i + 1) as usize] = l;
        }

        g.nvtxs = cnvtxs;
        g.nedges = l;
        g.ncon = 1;

        SetupGraph_tvwgt(&mut g);
        SetupGraph_label(&mut g);

        graph = Some(g);
    }

    // gk_free keys, map, mark: dropped (Vec Drop)

    graph
}

/// HiGHS: compress.c::PruneGraph (lines 152-222)
pub fn PruneGraph(
    ctrl: &mut ctrl_t,
    nvtxs: idx_t,
    xadj: &[idx_t],
    adjncy: &[idx_t],
    vwgt: Option<&[idx_t]>,
    iperm: &mut [idx_t],
    mut factor: real_t,
) -> Option<graph_t> {
    // IFSET(ctrl->dbglvl, METIS_DBG_INFO, ...) (compress.c 177-178) dropped

    let mut perm = vec![0 as idx_t; nvtxs as usize];

    // C: factor = factor*xadj[nvtxs]/nvtxs;
    // factor is real_t (f32). xadj[nvtxs] is idx_t (i32). nvtxs is idx_t (i32).
    // C: factor * xadj[nvtxs] -> f32 * i32 -> f32. Then / nvtxs -> f32 / i32 -> f32.
    // Result stored back into factor (f32).
    factor = (factor * (xadj[nvtxs as usize] as real_t)) / (nvtxs as real_t);

    let mut pnvtxs: idx_t = 0;
    let mut pnedges: idx_t = 0;
    #[allow(unused_assignments)]
    let mut nlarge: idx_t = 0;
    for i in 0..nvtxs {
        // C: xadj[i+1]-xadj[i] < factor
        // idx_t difference converts to real_t (f32) in mixed comparison.
        if ((xadj[(i + 1) as usize] - xadj[i as usize]) as real_t) < factor {
            perm[i as usize] = pnvtxs;
            iperm[pnvtxs as usize] = i;
            pnvtxs += 1;
            pnedges += xadj[(i + 1) as usize] - xadj[i as usize];
        } else {
            nlarge += 1;
            perm[i as usize] = nvtxs - nlarge;
            iperm[(nvtxs - nlarge) as usize] = i;
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_INFO, ...) (compress.c 177-178) dropped

    let mut graph: Option<graph_t> = None;

    if nlarge > 0 && nlarge < nvtxs {
        /* Prunning is possible, so go ahead and create the prunned graph */
        let mut g = CreateGraph();

        /* Allocate memory for the prunned graph */
        g.xadj = vec![0 as idx_t; (pnvtxs + 1) as usize];
        g.vwgt = vec![0 as idx_t; pnvtxs as usize];
        g.adjncy = vec![0 as idx_t; pnedges as usize];
        g.adjwgt = vec![1 as idx_t; pnedges as usize];

        g.xadj[0] = 0;
        let mut pnedges: idx_t = 0;
        let mut l: idx_t = 0;
        for i in 0..nvtxs {
            // C: xadj[i+1]-xadj[i] < factor
            // idx_t difference converts to real_t (f32) in mixed comparison.
            if ((xadj[(i + 1) as usize] - xadj[i as usize]) as real_t) < factor {
                // C: pvwgt[l] = (vwgt == NULL ? 1 : vwgt[i]);
                g.vwgt[l as usize] = match vwgt {
                    None => 1,
                    Some(v) => v[i as usize],
                };

                let mut j = xadj[i as usize];
                while j < xadj[(i + 1) as usize] {
                    let k = perm[adjncy[j as usize] as usize];
                    if k < pnvtxs {
                        g.adjncy[pnedges as usize] = k;
                        pnedges += 1;
                    }
                    j += 1;
                }
                l += 1;
                g.xadj[l as usize] = pnedges;
            }
        }

        g.nvtxs = pnvtxs;
        g.nedges = pnedges;
        g.ncon = 1;

        SetupGraph_tvwgt(&mut g);
        SetupGraph_label(&mut g);

        graph = Some(g);
    } else if nlarge > 0 && nlarge == nvtxs {
        // IFSET(ctrl->dbglvl, METIS_DBG_INFO, ...) (compress.c 213-214) dropped
        #[allow(unused_assignments)]
        {
            nlarge = 0;
        }
    }

    // gk_free perm: dropped (Vec Drop)

    graph
}

/// HiGHS: contig.c::FindSepInducedComponents (lines 21-83)
pub fn FindSepInducedComponents(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    cptr: &mut [idx_t],
    cind: &mut [idx_t],
) -> idx_t {
    let nvtxs = graph.nvtxs;
    let xadj = graph.xadj_eff();
    let adjncy = graph.adjncy_eff();
    let where_ = &graph.where_;

    let mut touched = vec![0 as idx_t; nvtxs as usize];

    for i in 0..graph.nbnd {
        touched[graph.bndind[i as usize] as usize] = 1;
    }

    // queue = cind; (alias)

    let mut nleft: idx_t = 0;
    for i in 0..nvtxs {
        if where_[i as usize] != 2 {
            nleft += 1;
        }
    }

    // C: for (i=0; i<nvtxs; i++) if (where[i] != 2) break;
    // i is left at nvtxs when every vertex is in the separator.
    let mut i: idx_t = 0;
    while i < nvtxs {
        if where_[i as usize] != 2 {
            break;
        }
        i += 1;
    }

    // C writes touched[nvtxs] (one past the ismalloc buffer, never read)
    // when i == nvtxs; that write is skipped here instead of panicking.
    if i < nvtxs {
        touched[i as usize] = 1;
    }
    cind[0] = i;
    let mut first: idx_t = 0;
    let mut last: idx_t = 1;
    cptr[0] = 0; /* This actually points to queue */
    let mut ncmps: idx_t = 0;

    while first != nleft {
        if first == last {
            /* Find another starting vertex */
            ncmps += 1;
            cptr[ncmps as usize] = first;
            // C: for (i=0; i<nvtxs; i++) if (!touched[i]) break;
            let mut ii: idx_t = 0;
            while ii < nvtxs {
                if touched[ii as usize] == 0 {
                    break;
                }
                ii += 1;
            }
            cind[last as usize] = ii;
            last += 1;
            touched[ii as usize] = 1;
        }

        let i = cind[first as usize];
        first += 1;
        let mut j = xadj[i as usize];
        while j < xadj[(i + 1) as usize] {
            let k = adjncy[j as usize];
            if touched[k as usize] == 0 {
                cind[last as usize] = k;
                last += 1;
                touched[k as usize] = 1;
            }
            j += 1;
        }
    }
    ncmps += 1;
    cptr[ncmps as usize] = first;

    // gk_free touched: dropped (Vec Drop)

    ncmps
}

/// HiGHS: separator.c::ConstructSeparator (lines 21-57)
pub fn ConstructSeparator(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // WCOREPUSH; // level-0 no-op, dropped

    let nvtxs = graph.nvtxs;
    let xadj = graph.xadj_eff();
    let nbnd = graph.nbnd;
    let bndind = &graph.bndind;

    // C: where = icopy(nvtxs, graph->where, iwspacemalloc(ctrl, nvtxs));
    // iwspacemalloc -> fresh Vec. icopy copies into it.
    let mut where_copy = vec![0 as idx_t; nvtxs as usize];
    icopy(nvtxs as usize, &graph.where_, &mut where_copy);

    /* Put the nodes in the boundary into the separator */
    for i in 0..nbnd {
        let j = bndind[i as usize];
        if xadj[(j + 1) as usize] - xadj[j as usize] > 0 {
            /* Ignore islands */
            where_copy[j as usize] = 2;
        }
    }

    FreeRData(graph);

    Allocate2WayNodePartitionMemory(ctrl, graph);
    icopy(nvtxs as usize, &where_copy, &mut graph.where_);

    // WCOREPOP; // level-0 no-op, dropped

    Compute2WayNodePartitionParams(ctrl, graph);

    FM_2WayNodeRefine2Sided(ctrl, graph, 1);
    FM_2WayNodeRefine1Sided(ctrl, graph, 4);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn compress_identical_pair() {
        // vertices 0 and 1 have identical closed neighbourhoods {0,1,2}
        let mut ctrl = ctrl_t {
            optype: crate::ipm::hipo::orderings::metis::METIS_OP_OMETIS,
            ..Default::default()
        };
        let xadj = [0, 2, 4, 6];
        let adjncy = [1, 2, 0, 2, 0, 1];
        let mut cptr = [0; 4];
        let mut cind = [0; 3];
        let g = CompressGraph(&mut ctrl, 3, &xadj, &adjncy, None, &mut cptr, &mut cind);
        assert!(g.is_none() || g.unwrap().nvtxs < 3);
    }

    #[test]
    fn compress_triangle_merges_all() {
        // K3: every closed neighbourhood is {0,1,2}, keys all 3 -> one vertex
        // (compress.c 46-77), 1 < 0.85*3 so the graph is built (compress.c 85).
        let mut ctrl = ctrl_t::default();
        let xadj = [0, 2, 4, 6];
        let adjncy = [1, 2, 0, 2, 0, 1];
        let mut cptr = [0; 4];
        let mut cind = [0; 3];
        let g = CompressGraph(&mut ctrl, 3, &xadj, &adjncy, None, &mut cptr, &mut cind)
            .expect("compressed");
        assert_eq!(g.nvtxs, 1);
        assert_eq!(g.nedges, 0);
        assert_eq!(g.xadj, vec![0, 0]);
        assert_eq!(g.vwgt, vec![3]);
        assert_eq!(g.tvwgt, vec![3]);
        assert_eq!(&cptr[..2], &[0, 3]);
    }

    #[test]
    fn compress_path_returns_none() {
        // Path 0-1-2: keys 1,3,3 but 1 and 2 differ in degree -> no merge,
        // cnvtxs == nvtxs is not < 0.85*nvtxs.
        let mut ctrl = ctrl_t::default();
        let xadj = [0, 1, 3, 4];
        let adjncy = [1, 0, 2, 1];
        let mut cptr = [0; 4];
        let mut cind = [0; 3];
        let vwgt = [1, 1, 1];
        let g = CompressGraph(
            &mut ctrl,
            3,
            &xadj,
            &adjncy,
            Some(&vwgt),
            &mut cptr,
            &mut cind,
        );
        assert!(g.is_none());
        assert_eq!(cptr, [0, 1, 2, 3]);
    }

    #[test]
    fn prune_star_center() {
        // Star with centre 0 and leaves 1..4: factor = 1.5*8/5 = 2.4, the
        // centre (degree 4) is pruned and placed last (compress.c 166-175).
        let mut ctrl = ctrl_t::default();
        let xadj = [0, 4, 5, 6, 7, 8];
        let adjncy = [1, 2, 3, 4, 0, 0, 0, 0];
        let vwgt = [5, 1, 2, 3, 4];
        let mut iperm = [0; 5];
        let g =
            PruneGraph(&mut ctrl, 5, &xadj, &adjncy, Some(&vwgt), &mut iperm, 1.5).expect("pruned");
        assert_eq!(iperm, [1, 2, 3, 4, 0]);
        assert_eq!(g.nvtxs, 4);
        assert_eq!(g.nedges, 0);
        assert_eq!(g.xadj, vec![0, 0, 0, 0, 0]);
        assert_eq!(g.vwgt, vec![1, 2, 3, 4]);
        assert_eq!(g.tvwgt, vec![10]);
    }

    #[test]
    fn prune_all_large_returns_none() {
        // K2 with factor 0.5: factor = 0.5*2/2 = 0.5, both degrees 1 >= 0.5,
        // nlarge == nvtxs -> pruning ignored (compress.c 211-215).
        let mut ctrl = ctrl_t::default();
        let xadj = [0, 1, 2];
        let adjncy = [1, 0];
        let mut iperm = [0; 2];
        let g = PruneGraph(&mut ctrl, 2, &xadj, &adjncy, None, &mut iperm, 0.5);
        assert!(g.is_none());
        assert_eq!(iperm, [1, 0]);
    }

    fn path_graph(where_: Vec<idx_t>, bndind: Vec<idx_t>) -> graph_t {
        let mut g = CreateGraph();
        g.nvtxs = 5;
        g.xadj = vec![0, 1, 3, 5, 7, 8];
        g.adjncy = vec![1, 0, 2, 1, 3, 2, 4, 3];
        g.where_ = where_;
        g.nbnd = bndind.len() as idx_t;
        g.bndind = bndind;
        g
    }

    #[test]
    fn sep_components_path() {
        // Path 0-1-2-3-4 with separator {2}: components {0,1} and {3,4}.
        let mut ctrl = ctrl_t::default();
        let mut g = path_graph(vec![0, 0, 2, 1, 1], vec![2]);
        let mut cptr = [0; 6];
        let mut cind = [0; 5];
        let ncmps = FindSepInducedComponents(&mut ctrl, &mut g, &mut cptr, &mut cind);
        assert_eq!(ncmps, 2);
        assert_eq!(&cptr[..3], &[0, 2, 4]);
        assert_eq!(&cind[..4], &[0, 1, 3, 4]);
    }

    #[test]
    fn sep_components_all_separator() {
        // Every vertex in the separator: nleft = 0, the start-vertex loop
        // leaves i = nvtxs (contig.c 46-49), one empty component.
        let mut ctrl = ctrl_t::default();
        let mut g = path_graph(vec![2; 5], vec![0, 1, 2, 3, 4]);
        let mut cptr = [9; 6];
        let mut cind = [9; 5];
        let ncmps = FindSepInducedComponents(&mut ctrl, &mut g, &mut cptr, &mut cind);
        assert_eq!(ncmps, 1);
        assert_eq!(&cptr[..2], &[0, 0]);
        assert_eq!(cind[0], 5);
    }

    fn node_partitioned_path7(ctrl: &mut ctrl_t, where_: &[idx_t]) -> graph_t {
        let mut g = CreateGraph();
        g.nvtxs = 7;
        g.ncon = 1;
        g.xadj = vec![0, 1, 3, 5, 7, 9, 11, 12];
        g.adjncy = vec![1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5];
        g.vwgt = vec![1; 7];
        SetupGraph_tvwgt(&mut g);
        Allocate2WayNodePartitionMemory(ctrl, &mut g);
        g.where_.copy_from_slice(where_);
        Compute2WayNodePartitionParams(ctrl, &mut g);
        g
    }

    #[test]
    fn node_balance_moves_separator() {
        // sfm.c::FM_2WayNodeBalance on path 0..6, where = [0,2,1,1,1,1,1]:
        // pwgts [1,5,1], mult 0.5, badmaxpwgt 3; the separator walks right
        // twice (vertices 1 and 2 move to side 0), vertex 3 is then skipped
        // because pwgts[0]+1 > badmaxpwgt (sfm.c 520-521).
        use super::super::sfm::FM_2WayNodeBalance;
        let mut ctrl = ctrl_t {
            ubfactors: vec![1.0],
            ..Default::default()
        };
        let mut g = node_partitioned_path7(&mut ctrl, &[0, 2, 1, 1, 1, 1, 1]);
        assert_eq!(g.pwgts, vec![1, 5, 1]);
        FM_2WayNodeBalance(&mut ctrl, &mut g);
        assert_eq!(g.where_, vec![0, 0, 0, 2, 1, 1, 1]);
        assert_eq!(g.pwgts, vec![3, 3, 1]);
        assert_eq!(g.mincut, 1);
        assert_eq!(g.nbnd, 1);
        assert_eq!(g.bndind[0], 3);
    }

    #[test]
    fn node_balance_already_balanced_returns() {
        // pwgts [3,3,1]: max 3 >= badmaxpwgt 3 but |3-3| < 3*7/7 -> early
        // return (sfm.c 480-481), nothing changes.
        use super::super::sfm::FM_2WayNodeBalance;
        let mut ctrl = ctrl_t {
            ubfactors: vec![1.0],
            ..Default::default()
        };
        let mut g = node_partitioned_path7(&mut ctrl, &[0, 0, 0, 2, 1, 1, 1]);
        let before = g.where_.clone();
        let rng = ctrl.rng_state;
        FM_2WayNodeBalance(&mut ctrl, &mut g);
        assert_eq!(g.where_, before);
        assert_eq!(g.pwgts, vec![3, 3, 1]);
        assert_eq!(ctrl.rng_state, rng);
    }
}
