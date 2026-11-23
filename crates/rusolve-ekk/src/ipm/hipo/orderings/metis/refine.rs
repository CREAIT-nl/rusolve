//! HiGHS source: extern/metis/libmetis/refine.c, extern/metis/libmetis/srefine.c, extern/metis/libmetis/bucketsort.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Edge-cut partition parameters (`refine.c`, used by the edge-based initial
//! separator), node-separator refinement driver and projection (`srefine.c`),
//! and the counting sort `BucketSortKeysInc`.

use super::{GraphId, GraphPool, ctrl_t, graph_t, idx_t};

/// HiGHS: refine.c::Allocate2WayPartitionMemory (lines 17-30)
pub fn Allocate2WayPartitionMemory(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // HiGHS: refine.c lines 19-22: idx_t nvtxs, ncon; nvtxs = graph->nvtxs; ncon = graph->ncon;
    let nvtxs = graph.nvtxs;
    let ncon = graph.ncon;

    // C lines 24-29: graph->pwgts = imalloc(2*ncon); graph->where = imalloc(nvtxs); ...
    graph.pwgts = super::gklib::imalloc((2 * ncon) as usize);
    graph.where_ = super::gklib::imalloc(nvtxs as usize);
    graph.bndptr = super::gklib::imalloc(nvtxs as usize);
    graph.bndind = super::gklib::imalloc(nvtxs as usize);
    graph.id = super::gklib::imalloc(nvtxs as usize);
    graph.ed = super::gklib::imalloc(nvtxs as usize);
}

/// HiGHS: refine.c::Compute2WayPartitionParams (lines 36-99)
pub fn Compute2WayPartitionParams(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // HiGHS: refine.c lines 38-41: idx_t i, j, nvtxs, ncon, nbnd, mincut, istart, iend, tid, ted, me;
    // C lines 43-48:
    let nvtxs = graph.nvtxs;
    let ncon = graph.ncon;
    // C line 45: xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
    // C line 47: adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
    // vwgt = graph->vwgt; adjwgt = graph->adjwgt;
    // where = graph->where; id = graph->id; ed = graph->ed;

    // C line 54: pwgts = iset(2*ncon, 0, graph->pwgts);
    super::gklib::iset((2 * ncon) as usize, 0, &mut graph.pwgts);
    // C line 55: bndptr = iset(nvtxs, -1, graph->bndptr);
    super::gklib::iset(nvtxs as usize, -1, &mut graph.bndptr);
    // C line 56: bndind = graph->bndind;

    // HiGHS: refine.c lines 58-70: Compute pwgts
    if ncon == 1 {
        // C lines 60-62: for (i=0; i<nvtxs; i++) { pwgts[where[i]] += vwgt[i]; }
        for i in 0..nvtxs {
            let me = graph.where_[i as usize];
            graph.pwgts[me as usize] += graph.vwgt[i as usize];
        }
    } else {
        // C lines 65-69: for (i=0; i<nvtxs; i++) { me = where[i]; for (j=0; j<ncon; j++) pwgts[me*ncon+j] += vwgt[i*ncon+j]; }
        for i in 0..nvtxs {
            let me = graph.where_[i as usize];
            for j in 0..ncon {
                graph.pwgts[(me * ncon + j) as usize] += graph.vwgt[(i * ncon + j) as usize];
            }
        }
    }

    // HiGHS: refine.c lines 73-94: Compute the required info for refinement
    // C line 74: for (nbnd=0, mincut=0, i=0; i<nvtxs; i++)
    let mut nbnd: idx_t = 0;
    let mut mincut: idx_t = 0;
    for i in 0..nvtxs {
        // C line 75: istart = xadj[i];
        let istart = graph.xadj_eff()[i as usize];
        // C line 76: iend = xadj[i+1];
        let iend = graph.xadj_eff()[(i + 1) as usize];

        // C line 78: me = where[i];
        let me = graph.where_[i as usize];
        // C line 79: tid = ted = 0;
        let mut tid: idx_t = 0;
        let mut ted: idx_t = 0;

        // C lines 81-86: for (j=istart; j<iend; j++) { if (me == where[adjncy[j]]) tid += adjwgt[j]; else ted += adjwgt[j]; }
        for j in istart..iend {
            if me == graph.where_[graph.adjncy_eff()[j as usize] as usize] {
                tid += graph.adjwgt[j as usize];
            } else {
                ted += graph.adjwgt[j as usize];
            }
        }
        // C lines 87-88: id[i] = tid; ed[i] = ted;
        graph.id[i as usize] = tid;
        graph.ed[i as usize] = ted;

        // C lines 90-93: if (ted > 0 || istart == iend) { BNDInsert(nbnd, bndind, bndptr, i); mincut += ted; }
        if ted > 0 || istart == iend {
            super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, i);
            mincut += ted;
        }
    }

    // HiGHS: refine.c lines 96-97: graph->mincut = mincut/2; graph->nbnd = nbnd;
    graph.mincut = mincut / 2;
    graph.nbnd = nbnd;
}

/// HiGHS: srefine.c::Refine2WayNode (lines 23-53)
/// Calls `FM_2WayNodeBalance` and the `FM_2WayNodeRefine*` routines.
/// Walks `graph` up the `finer` chain to `orggraph`.
pub fn Refine2WayNode(ctrl: &mut ctrl_t, pool: &mut GraphPool, orggraph: GraphId, graph: GraphId) {
    // HiGHS: srefine.c line 27: if (graph == orggraph)
    if graph == orggraph {
        // HiGHS: srefine.c line 28: Compute2WayNodePartitionParams(ctrl, graph);
        {
            let g = pool.get_mut(graph);
            Compute2WayNodePartitionParams(ctrl, g);
        }
    } else {
        // HiGHS: srefine.c lines 31-50: do { ... } while (graph != orggraph);
        let mut current = graph;
        loop {
            // HiGHS: srefine.c line 32: graph = graph->finer;
            current = pool
                .get(current)
                .finer
                .expect("Refine2WayNode: graph->finer is NULL");

            // HiGHS: srefine.c line 34: Project2WayNodePartition(ctrl, graph);
            Project2WayNodePartition(ctrl, pool, current);

            // HiGHS: srefine.c line 36: FM_2WayNodeBalance(ctrl, graph);
            {
                let g = pool.get_mut(current);
                super::sfm::FM_2WayNodeBalance(ctrl, g);
            }

            // HiGHS: srefine.c lines 39-48: switch (ctrl->rtype)
            match ctrl.rtype {
                super::METIS_RTYPE_SEP2SIDED => {
                    // HiGHS: srefine.c line 41: FM_2WayNodeRefine2Sided(ctrl, graph, ctrl->niter);
                    let niter = ctrl.niter;
                    let g = pool.get_mut(current);
                    super::sfm::FM_2WayNodeRefine2Sided(ctrl, g, niter);
                }
                super::METIS_RTYPE_SEP1SIDED => {
                    // HiGHS: srefine.c line 44: FM_2WayNodeRefine1Sided(ctrl, graph, ctrl->niter);
                    let niter = ctrl.niter;
                    let g = pool.get_mut(current);
                    super::sfm::FM_2WayNodeRefine1Sided(ctrl, g, niter);
                }
                _ => {
                    // HiGHS: srefine.c line 47: GK_ERREXIT("Unknown rtype of %d\n", ctrl->rtype);
                    super::gklib::gk_errexit(&format!("Unknown rtype of {}\n", ctrl.rtype));
                }
            }

            // HiGHS: srefine.c line 50: } while (graph != orggraph);
            if current == orggraph {
                break;
            }
        }
    }
}

/// HiGHS: srefine.c::Allocate2WayNodePartitionMemory (lines 59-70)
pub fn Allocate2WayNodePartitionMemory(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // HiGHS: srefine.c lines 61-63: idx_t nvtxs; nvtxs = graph->nvtxs;
    let nvtxs = graph.nvtxs;

    // C line 65: graph->pwgts = imalloc(3);
    graph.pwgts = super::gklib::imalloc(3);
    // C line 66: graph->where = imalloc(nvtxs);
    graph.where_ = super::gklib::imalloc(nvtxs as usize);
    // C line 67: graph->bndptr = imalloc(nvtxs);
    graph.bndptr = super::gklib::imalloc(nvtxs as usize);
    // C line 68: graph->bndind = imalloc(nvtxs);
    graph.bndind = super::gklib::imalloc(nvtxs as usize);
    // C line 69: graph->nrinfo = (nrinfo_t *)malloc(nvtxs*sizeof(nrinfo_t));
    graph.nrinfo = vec![super::nrinfo_t::default(); nvtxs as usize];
}

/// HiGHS: srefine.c::Compute2WayNodePartitionParams (lines 76-122)
pub fn Compute2WayNodePartitionParams(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // HiGHS: srefine.c lines 78-83: idx_t i, j, nvtxs, nbnd; ...
    // C lines 85-88:
    let nvtxs = graph.nvtxs;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
    // vwgt = graph->vwgt;
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;

    // C lines 90-94:
    // where = graph->where;
    // rinfo = graph->nrinfo;
    // pwgts = iset(3, 0, graph->pwgts);
    super::gklib::iset(3, 0, &mut graph.pwgts);
    // bndind = graph->bndind;
    // bndptr = iset(nvtxs, -1, graph->bndptr);
    super::gklib::iset(nvtxs as usize, -1, &mut graph.bndptr);

    // HiGHS: srefine.c lines 97-117: Compute now the separator external degrees
    // C line 100: nbnd = 0;
    let mut nbnd: idx_t = 0;
    // C line 101: for (i=0; i<nvtxs; i++)
    for i in 0..nvtxs {
        // C line 102: me = where[i];
        let me = graph.where_[i as usize];
        // C line 103: pwgts[me] += vwgt[i];
        graph.pwgts[me as usize] += graph.vwgt[i as usize];

        // C line 105: if (me == 2) { ... }
        if me == 2 {
            // C line 106: BNDInsert(nbnd, bndind, bndptr, i);
            super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, i);

            // C line 108: edegrees = rinfo[i].edegrees;
            // C line 109: edegrees[0] = edegrees[1] = 0;
            graph.nrinfo[i as usize].edegrees[0] = 0;
            graph.nrinfo[i as usize].edegrees[1] = 0;

            // C lines 111-115: for (j=xadj[i]; j<xadj[i+1]; j++) { other = where[adjncy[j]]; if (other != 2) edegrees[other] += vwgt[adjncy[j]]; }
            for j in graph.xadj_eff()[i as usize]..graph.xadj_eff()[(i + 1) as usize] {
                let other = graph.where_[graph.adjncy_eff()[j as usize] as usize];
                if other != 2 {
                    graph.nrinfo[i as usize].edegrees[other as usize] +=
                        graph.vwgt[graph.adjncy_eff()[j as usize] as usize];
                }
            }
        }
    }

    // HiGHS: srefine.c lines 120-121: graph->mincut = pwgts[2]; graph->nbnd = nbnd;
    graph.mincut = graph.pwgts[2];
    graph.nbnd = nbnd;
}

/// HiGHS: srefine.c::Project2WayNodePartition (lines 128-152)
/// Projects from `graph->coarser` and frees it.
pub fn Project2WayNodePartition(ctrl: &mut ctrl_t, pool: &mut GraphPool, graph: GraphId) {
    // HiGHS: srefine.c lines 130-131: idx_t i, j, nvtxs; idx_t *cmap, *where, *cwhere; graph_t *cgraph;
    // C line 134: cgraph = graph->coarser;
    // C line 135: cwhere = cgraph->where;
    // C line 137: nvtxs = graph->nvtxs;
    // C line 138: cmap = graph->cmap;

    // C line 140: Allocate2WayNodePartitionMemory(ctrl, graph);
    {
        let g = pool.get_mut(graph);
        Allocate2WayNodePartitionMemory(ctrl, g);
    }

    // C line 141: where = graph->where;
    // C lines 143-146: for (i=0; i<nvtxs; i++) { where[i] = cwhere[cmap[i]]; }
    {
        let cgraph_id = pool
            .get(graph)
            .coarser
            .expect("Project2WayNodePartition: graph->coarser is NULL");
        let (g, cg) = pool.get2_mut(graph, cgraph_id);
        let nvtxs = g.nvtxs;
        for i in 0..nvtxs {
            g.where_[i as usize] = cg.where_[g.cmap[i as usize] as usize];
        }
    }

    // C line 148: FreeGraph(&graph->coarser);
    // C line 149: graph->coarser = NULL;
    {
        let mut coarser = pool.get_mut(graph).coarser;
        super::graph::FreeGraph(pool, &mut coarser);
        pool.get_mut(graph).coarser = None;
    }

    // C line 151: Compute2WayNodePartitionParams(ctrl, graph);
    {
        let g = pool.get_mut(graph);
        Compute2WayNodePartitionParams(ctrl, g);
    }
}

/// HiGHS: bucketsort.c::BucketSortKeysInc (lines 23-43)
pub fn BucketSortKeysInc(
    ctrl: &mut ctrl_t,
    n: idx_t,
    max: idx_t,
    keys: &[idx_t],
    tperm: &[idx_t],
    perm: &mut [idx_t],
) {
    // HiGHS: bucketsort.c line 29: WCOREPUSH;
    super::wspace::wspacepush(ctrl);

    // HiGHS: bucketsort.c line 31: counts = iset(max+2, 0, iwspacemalloc(ctrl, max+2));
    let mut counts = super::wspace::iwspacemalloc(ctrl, max + 2);
    super::gklib::iset((max + 2) as usize, 0, &mut counts);

    // HiGHS: bucketsort.c lines 33-34: for (i=0; i<n; i++) counts[keys[i]]++;
    for i in 0..n {
        counts[keys[i as usize] as usize] += 1;
    }

    // HiGHS: bucketsort.c line 35: MAKECSR(i, max+1, counts);
    // GKlib/gk_macros.h lines 41-46: MAKECSR(i, n, a) expands to
    //   for (i=1; i<n; i++) a[i] += a[i-1];
    //   for (i=n; i>0; i--) a[i] = a[i-1];
    //   a[0] = 0;
    {
        let n = max + 1;
        let mut i: idx_t = 1;
        while i < n {
            counts[i as usize] += counts[(i - 1) as usize];
            i += 1;
        }
        i = n;
        while i > 0 {
            counts[i as usize] = counts[(i - 1) as usize];
            i -= 1;
        }
        counts[0] = 0;
    }

    // HiGHS: bucketsort.c lines 37-40: for (ii=0; ii<n; ii++) { i = tperm[ii]; perm[counts[keys[i]]++] = i; }
    for ii in 0..n {
        let i = tperm[ii as usize];
        perm[counts[keys[i as usize] as usize] as usize] = i;
        counts[keys[i as usize] as usize] += 1;
    }

    // HiGHS: bucketsort.c line 42: WCOREPOP;
    super::wspace::wspacepop(ctrl);
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn bucket_sort_stable_increasing() {
        let mut ctrl = ctrl_t::default();
        let keys = [2, 0, 1, 0];
        let tperm = [0, 1, 2, 3];
        let mut perm = [0; 4];
        BucketSortKeysInc(&mut ctrl, 4, 2, &keys, &tperm, &mut perm);
        assert_eq!(perm, [1, 3, 2, 0]);
    }

    #[test]
    fn bucket_sort_makecsr_prefix_sum() {
        // Counts [1, 2, 1] need the MAKECSR prefix sum to place every key.
        let mut ctrl = ctrl_t::default();
        let keys = [1, 1, 0, 2];
        let tperm = [0, 1, 2, 3];
        let mut perm = [-1; 4];
        BucketSortKeysInc(&mut ctrl, 4, 2, &keys, &tperm, &mut perm);
        assert_eq!(perm, [2, 0, 1, 3]);
    }
}
