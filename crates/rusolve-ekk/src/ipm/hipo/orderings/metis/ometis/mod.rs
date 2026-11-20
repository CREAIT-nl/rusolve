//! HiGHS source: extern/metis/libmetis/ometis.c (lines 1-306)
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! `Highs_METIS_NodeND` and the nested-dissection recursion. Graphs are
//! created and freed in the [`super::GraphPool`] of the call; the recursion
//! order (left subgraph first, then right) fixes the `rng_state` sequence.

pub mod part2;

pub use part2::*;

use super::{GraphId, GraphPool, ctrl_t, graph_t, idx_t};

/// HiGHS: ometis.c::Highs_METIS_NodeND (lines 43-140)
/// `vwgt` is `None` for the C `NULL` (HiPO passes `nullptr`); `options`
/// is the `METIS_NOPTIONS` array. Returns `METIS_OK` / `METIS_ERROR_INPUT`.
pub fn Highs_METIS_NodeND(
    nvtxs: &mut idx_t,
    xadj: &[idx_t],
    adjncy: &[idx_t],
    vwgt: Option<&[idx_t]>,
    options: Option<&[idx_t]>,
    perm: &mut [idx_t],
    iperm: &mut [idx_t],
) -> i32 {
    // ometis.c lines 46-49: locals
    let mut nnvtxs: idx_t = 0;
    let mut pool = GraphPool::new();
    // `graph` is the C `graph_t *graph = NULL`; the pool id replaces the pointer.
    let mut graph: Option<GraphId> = None;
    let mut piperm: Vec<idx_t> = Vec::new();
    let mut cptr: Vec<idx_t> = Vec::new();
    let mut cind: Vec<idx_t> = Vec::new();

    // ometis.c line 52: ctrl = SetupCtrl(METIS_OP_OMETIS, options, 1, 3, NULL, NULL);
    // ometis.c lines 53-55: if (!ctrl) return METIS_ERROR_INPUT;
    let mut ctrl: Box<ctrl_t> =
        match super::options::SetupCtrl(super::METIS_OP_OMETIS, options, 1, 3, None, None) {
            Some(c) => c,
            None => return super::METIS_ERROR_INPUT,
        };

    // ometis.c lines 57-71: prune the dense columns
    if ctrl.pfactor > 0.0 {
        // ometis.c line 59: piperm = imalloc(*nvtxs);
        piperm = vec![0; *nvtxs as usize];

        // ometis.c line 61: graph = PruneGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, piperm, ctrl->pfactor);
        let pfactor = ctrl.pfactor;
        match super::compress::PruneGraph(
            &mut ctrl,
            *nvtxs,
            xadj,
            adjncy,
            vwgt,
            &mut piperm,
            pfactor,
        ) {
            None => {
                // ometis.c lines 63-65: gk_free(&piperm) (Drop); ctrl->pfactor = 0.0;
                ctrl.pfactor = 0.0;
            }
            Some(gg) => {
                // ometis.c lines 68-69: nnvtxs = graph->nvtxs; ctrl->compress = 0;
                nnvtxs = gg.nvtxs;
                graph = Some(pool.alloc(gg));
                ctrl.compress = 0;
            }
        }
    }

    // ometis.c lines 73-93: compress the graph
    if ctrl.compress != 0 {
        // ometis.c line 76: cptr = imalloc(*nvtxs+1);
        cptr = vec![0; (*nvtxs + 1) as usize];
        // ometis.c line 77: cind = imalloc(*nvtxs);
        cind = vec![0; *nvtxs as usize];

        // ometis.c line 79: graph = CompressGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, cptr, cind);
        match super::compress::CompressGraph(
            &mut ctrl, *nvtxs, xadj, adjncy, vwgt, &mut cptr, &mut cind,
        ) {
            None => {
                // ometis.c lines 82-84: gk_free(&cptr), gk_free(&cind) (Drop); ctrl->compress = 0;
                ctrl.compress = 0;
            }
            Some(gg) => {
                // ometis.c lines 87-90
                nnvtxs = gg.nvtxs;
                graph = Some(pool.alloc(gg));
                // ometis.c line 88: ctrl->cfactor = 1.0*(*nvtxs)/nnvtxs;
                // double expression, stored into real_t (f32).
                ctrl.cfactor = ((1.0 * (*nvtxs as f64)) / (nnvtxs as f64)) as f32;
                // ometis.c line 89: if (ctrl->cfactor > 1.5 && ctrl->nseps == 1)
                if (ctrl.cfactor as f64) > 1.5 && ctrl.nseps == 1 {
                    // ometis.c line 90
                    ctrl.nseps = 2;
                }
                // ometis.c line 91: //ctrl->nseps = (idx_t)(ctrl->cfactor*ctrl->nseps); dropped (comment)
            }
        }
    }

    // ometis.c lines 95-97: if no prunning and no compression, setup the graph in the normal way.
    if ctrl.pfactor == 0.0 && ctrl.compress == 0 {
        let gg = super::graph::SetupGraph(&mut ctrl, *nvtxs, 1, xadj, adjncy, vwgt, None, None);
        graph = Some(pool.alloc(gg));
    }

    // Every path above leaves `graph` set (C dereferences `graph->nvtxs` below).
    let graph_id = graph.expect("ometis.c: graph is NULL after setup");

    // ometis.c line 101: AllocateWorkSpace(ctrl, graph);
    super::wspace::AllocateWorkSpace(&mut ctrl, pool.get(graph_id));

    // ometis.c lines 104-107: do the nested dissection ordering
    if ctrl.ccorder != 0 {
        // ometis.c line 105: MlevelNestedDissectionCC(ctrl, graph, iperm, graph->nvtxs);
        let gnvtxs = pool.get(graph_id).nvtxs;
        MlevelNestedDissectionCC(&mut ctrl, &mut pool, graph_id, iperm, gnvtxs);
    } else {
        // ometis.c line 107: MlevelNestedDissection(ctrl, graph, iperm, graph->nvtxs);
        let gnvtxs = pool.get(graph_id).nvtxs;
        MlevelNestedDissection(&mut ctrl, &mut pool, graph_id, iperm, gnvtxs);
    }

    // ometis.c lines 110-118: Order any prunned vertices
    if ctrl.pfactor > 0.0 {
        // ometis.c line 111: icopy(nnvtxs, iperm, perm);  /* Use perm as an auxiliary array */
        super::gklib::icopy(nnvtxs as usize, iperm, perm);
        // ometis.c lines 112-113
        for i in 0..nnvtxs {
            iperm[piperm[i as usize] as usize] = perm[i as usize];
        }
        // ometis.c lines 114-115
        for i in nnvtxs..*nvtxs {
            iperm[piperm[i as usize] as usize] = i;
        }
        // ometis.c line 117: gk_free((void **)&piperm); (Drop)
    } else if ctrl.compress != 0 {
        // ometis.c lines 119-131: Uncompress the ordering
        // ometis.c lines 121-122: construct perm from iperm
        for i in 0..nnvtxs {
            perm[iperm[i as usize] as usize] = i;
        }
        // ometis.c lines 123-127
        let mut l: idx_t = 0;
        for ii in 0..nnvtxs {
            let i = perm[ii as usize];
            let mut j = cptr[i as usize];
            while j < cptr[(i + 1) as usize] {
                iperm[cind[j as usize] as usize] = l;
                l += 1;
                j += 1;
            }
        }
        // ometis.c lines 129-130: gk_free(&cptr); gk_free(&cind); (Drop)
    }

    // ometis.c lines 133-134
    for i in 0..*nvtxs {
        perm[iperm[i as usize] as usize] = i;
    }

    // ometis.c line 137: FreeCtrl(&ctrl);
    let mut ctrl_opt = Some(ctrl);
    super::options::FreeCtrl(&mut ctrl_opt);

    // ometis.c line 139: return METIS_OK;
    super::METIS_OK
}

/// HiGHS: ometis.c::MlevelNestedDissection (lines 151-193)
pub fn MlevelNestedDissection(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    order: &mut [idx_t],
    lastvtx: idx_t,
) {
    // ometis.c line 158: nvtxs = graph->nvtxs; (only read by the dropped IFSET debug print)

    // ometis.c line 160: MlevelNodeBisectionMultiple(ctrl, graph);
    MlevelNodeBisectionMultiple(ctrl, pool, graph);

    // ometis.c lines 162-164: IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO, ...) dropped (debug print)

    // ometis.c lines 168-172: Order the nodes in the separator
    let nbnd = pool.get(graph).nbnd;
    let mut lastvtx = lastvtx;
    for i in 0..nbnd {
        // ometis.c line 172: order[label[bndind[i]]] = --lastvtx;
        lastvtx -= 1;
        let bi = pool.get(graph).bndind[i as usize];
        let lb = pool.get(graph).label[bi as usize];
        order[lb as usize] = lastvtx;
    }

    // ometis.c line 174: SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);
    let mut lgraph_opt: Option<GraphId> = None;
    let mut rgraph_opt: Option<GraphId> = None;
    SplitGraphOrder(ctrl, pool, graph, &mut lgraph_opt, &mut rgraph_opt);
    let lgraph = lgraph_opt.expect("SplitGraphOrder: lgraph NULL");
    let rgraph = rgraph_opt.expect("SplitGraphOrder: rgraph NULL");

    // ometis.c line 177: FreeGraph(&graph);
    let mut graph_opt = Some(graph);
    super::graph::FreeGraph(pool, &mut graph_opt);

    // ometis.c lines 181-186: Recurse on lgraph first, as its lastvtx depends on rgraph->nvtxs
    let rgraph_nvtxs = pool.get(rgraph).nvtxs;
    if pool.get(lgraph).nvtxs > super::MMDSWITCH && pool.get(lgraph).nedges > 0 {
        MlevelNestedDissection(ctrl, pool, lgraph, order, lastvtx - rgraph_nvtxs);
    } else {
        // ometis.c line 184: MMDOrder(ctrl, lgraph, order, lastvtx-rgraph->nvtxs);
        MMDOrder(ctrl, pool.get_mut(lgraph), order, lastvtx - rgraph_nvtxs);
        // ometis.c line 185: FreeGraph(&lgraph);
        let mut lgraph_opt = Some(lgraph);
        super::graph::FreeGraph(pool, &mut lgraph_opt);
    }

    // ometis.c lines 187-192
    if pool.get(rgraph).nvtxs > super::MMDSWITCH && pool.get(rgraph).nedges > 0 {
        MlevelNestedDissection(ctrl, pool, rgraph, order, lastvtx);
    } else {
        // ometis.c line 190: MMDOrder(ctrl, rgraph, order, lastvtx);
        MMDOrder(ctrl, pool.get_mut(rgraph), order, lastvtx);
        // ometis.c line 191: FreeGraph(&rgraph);
        let mut rgraph_opt = Some(rgraph);
        super::graph::FreeGraph(pool, &mut rgraph_opt);
    }
}

/// HiGHS: ometis.c::MlevelNestedDissectionCC (lines 204-261)
pub fn MlevelNestedDissectionCC(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    order: &mut [idx_t],
    lastvtx: idx_t,
) {
    // ometis.c line 212: nvtxs = graph->nvtxs;
    let nvtxs = pool.get(graph).nvtxs;

    // ometis.c line 214: MlevelNodeBisectionMultiple(ctrl, graph);
    MlevelNodeBisectionMultiple(ctrl, pool, graph);

    // ometis.c lines 216-218: IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO, ...) dropped (debug print)

    // ometis.c lines 220-225: Order the nodes in the separator
    let nbnd = pool.get(graph).nbnd;
    let mut lastvtx = lastvtx;
    for i in 0..nbnd {
        lastvtx -= 1;
        let bi = pool.get(graph).bndind[i as usize];
        let lb = pool.get(graph).label[bi as usize];
        order[lb as usize] = lastvtx;
    }

    // ometis.c line 227: WCOREPUSH;
    super::wspace::wspacepush(ctrl);

    // ometis.c lines 228-229
    let mut cptr = super::wspace::iwspacemalloc(ctrl, nvtxs + 1);
    let mut cind = super::wspace::iwspacemalloc(ctrl, nvtxs);

    // ometis.c line 230: ncmps = FindSepInducedComponents(ctrl, graph, cptr, cind);
    let ncmps =
        super::compress::FindSepInducedComponents(ctrl, pool.get_mut(graph), &mut cptr, &mut cind);

    // ometis.c lines 232-235: IFSET(ctrl->dbglvl&METIS_DBG_INFO, ...) dropped (debug print)

    // ometis.c line 237: sgraphs = SplitGraphOrderCC(ctrl, graph, ncmps, cptr, cind);
    let sgraphs: Vec<GraphId> = SplitGraphOrderCC(ctrl, pool, graph, ncmps, &cptr, &cind);

    // ometis.c line 239: WCOREPOP;
    super::wspace::wspacepop(ctrl);

    // ometis.c lines 241-242: FreeGraph(&graph);
    let mut graph_opt = Some(graph);
    super::graph::FreeGraph(pool, &mut graph_opt);

    // ometis.c lines 245-258: Go and process the subgraphs
    let mut rnvtxs: idx_t = 0;
    for i in 0..ncmps {
        let sg = sgraphs[i as usize];
        // ometis.c lines 246-248: snvtxs = sgraphs[i]->nvtxs;
        let snvtxs = pool.get(sg).nvtxs;

        // ometis.c lines 250-256
        if pool.get(sg).nvtxs > super::MMDSWITCH && pool.get(sg).nedges > 0 {
            MlevelNestedDissectionCC(ctrl, pool, sg, order, lastvtx - rnvtxs);
        } else {
            MMDOrder(ctrl, pool.get_mut(sg), order, lastvtx - rnvtxs);
            let mut sg_opt = Some(sg);
            super::graph::FreeGraph(pool, &mut sg_opt);
        }
        // ometis.c line 257
        rnvtxs += snvtxs;
    }

    // ometis.c line 260: gk_free((void **)&sgraphs); (Drop)
}

/// HiGHS: ometis.c::MlevelNodeBisectionMultiple (lines 268-306)
pub fn MlevelNodeBisectionMultiple(ctrl: &mut ctrl_t, pool: &mut GraphPool, graph: GraphId) {
    // ometis.c lines 273-277: if the graph is small, just find a single vertex separator
    if ctrl.nseps == 1 || pool.get(graph).nvtxs < (if ctrl.compress != 0 { 1000 } else { 2000 }) {
        MlevelNodeBisectionL2(ctrl, pool, graph, super::LARGENIPARTS);
        return;
    }

    // ometis.c line 279: WCOREPUSH;
    super::wspace::wspacepush(ctrl);

    // ometis.c line 281: bestwhere = iwspacemalloc(ctrl, graph->nvtxs);
    let mut bestwhere = super::wspace::iwspacemalloc(ctrl, pool.get(graph).nvtxs);

    // ometis.c line 283: mincut = graph->tvwgt[0];
    let mut mincut = pool.get(graph).tvwgt[0];

    // ometis.c lines 284-298
    let mut i: idx_t = 0;
    while i < ctrl.nseps {
        // ometis.c line 285: MlevelNodeBisectionL2(ctrl, graph, LARGENIPARTS);
        MlevelNodeBisectionL2(ctrl, pool, graph, super::LARGENIPARTS);

        // ometis.c lines 287-291
        if i == 0 || pool.get(graph).mincut < mincut {
            mincut = pool.get(graph).mincut;
            if i < ctrl.nseps - 1 {
                // ometis.c line 290: icopy(graph->nvtxs, graph->where, bestwhere);
                let g = pool.get(graph);
                super::gklib::icopy(g.nvtxs as usize, &g.where_, &mut bestwhere);
            }
        }

        // ometis.c lines 293-294
        if mincut == 0 {
            break;
        }

        // ometis.c lines 296-297: if (i < ctrl->nseps-1) FreeRData(graph);
        if i < ctrl.nseps - 1 {
            super::graph::FreeRData(pool.get_mut(graph));
        }

        i += 1;
    }

    // ometis.c lines 300-303
    if mincut != pool.get(graph).mincut {
        // ometis.c line 301: icopy(graph->nvtxs, bestwhere, graph->where);
        let g = pool.get_mut(graph);
        super::gklib::icopy(g.nvtxs as usize, &bestwhere, &mut g.where_);
        // ometis.c line 302: Compute2WayNodePartitionParams(ctrl, graph);
        super::refine::Compute2WayNodePartitionParams(ctrl, g);
    }

    // ometis.c line 305: WCOREPOP;
    super::wspace::wspacepop(ctrl);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn node_nd_small_path() {
        // path on 4 vertices; any valid permutation is accepted here, the
        // exact one is checked against the oracle in t_w4_orderings.
        let mut n: idx_t = 4;
        let xadj = [0, 1, 3, 5, 6];
        let adjncy = [1, 0, 2, 1, 3, 2];
        let mut options = vec![-1 as idx_t; crate::ipm::hipo::orderings::metis::METIS_NOPTIONS];
        options[crate::ipm::hipo::orderings::metis::METIS_OPTION_SEED] = 42;
        options[crate::ipm::hipo::orderings::metis::METIS_OPTION_NO2HOP] = 1;
        let (mut perm, mut iperm) = (vec![0; 4], vec![0; 4]);
        let st = Highs_METIS_NodeND(
            &mut n,
            &xadj,
            &adjncy,
            None,
            Some(&options),
            &mut perm,
            &mut iperm,
        );
        assert_eq!(st, crate::ipm::hipo::orderings::metis::METIS_OK);
        let mut p = perm.clone();
        p.sort();
        assert_eq!(p, vec![0, 1, 2, 3]);
        for i in 0..4 {
            assert_eq!(perm[iperm[i] as usize] as usize, i);
        }
    }
}
