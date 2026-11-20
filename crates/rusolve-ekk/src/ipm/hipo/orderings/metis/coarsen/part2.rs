//! HiGHS source: extern/metis/libmetis/coarsen.c (lines 416-895)
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! 2-hop matchings and coarse graph contraction (`CreateCoarseGraph` with its
//! `HTLENGTH` hash table and the direct-access path for large degrees).

use super::super::gklib::{
    iaxpy, icopy, ikv_t, ikvsorti, imalloc, irandInRange, iset, isortd, isum, rmalloc,
};
use super::super::graph::CreateGraph;
use super::super::wspace::{ikvwspacemalloc, iwspacemalloc, wspacepop, wspacepush};
use super::super::{
    GraphId, GraphPool, HTLENGTH, IDX_MAX, UNMATCHED, UNMATCHEDFOR2HOP, ctrl_t, graph_t, idx_t,
    real_t,
};
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: coarsen.c::Match_2Hop (lines 422-434)
pub fn Match_2Hop(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    perm: &mut [idx_t],
    match_: &mut [idx_t],
    cnvtxs: idx_t,
    nunmatched: usize,
) -> idx_t {
    let mut nunmatched = nunmatched;
    let mut cnvtxs = Match_2HopAny(ctrl, graph, perm, match_, cnvtxs, &mut nunmatched, 2);
    cnvtxs = Match_2HopAll(ctrl, graph, perm, match_, cnvtxs, &mut nunmatched, 64);
    if (nunmatched as f64) > 1.5 * UNMATCHEDFOR2HOP * (graph.nvtxs as f64) {
        cnvtxs = Match_2HopAny(ctrl, graph, perm, match_, cnvtxs, &mut nunmatched, 3);
    }
    if (nunmatched as f64) > 2.0 * UNMATCHEDFOR2HOP * (graph.nvtxs as f64) {
        let maxdegree = graph.nvtxs as usize;
        cnvtxs = Match_2HopAny(
            ctrl,
            graph,
            perm,
            match_,
            cnvtxs,
            &mut nunmatched,
            maxdegree,
        );
    }
    cnvtxs
}

/// HiGHS: coarsen.c::Match_2HopAny (lines 444-506)
pub fn Match_2HopAny(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    perm: &mut [idx_t],
    match_: &mut [idx_t],
    cnvtxs: idx_t,
    r_nunmatched: &mut usize,
    maxdegree: usize,
) -> idx_t {
    let mut cnvtxs = cnvtxs;
    let nvtxs = graph.nvtxs;
    let xadj: &[idx_t] = match &graph.xadj_orig {
        Some(x) => x.as_slice(),
        None => graph.xadj.as_slice(),
    };
    let adjncy: &[idx_t] = match &graph.adjncy_orig {
        Some(a) => a.as_slice(),
        None => graph.adjncy.as_slice(),
    };
    let cmap = &mut graph.cmap;

    let mut nunmatched = *r_nunmatched;

    // create the inverted index
    wspacepush(ctrl);
    let mut colptr = iwspacemalloc(ctrl, nvtxs + 1);
    iset(nvtxs as usize, 0, &mut colptr);
    for i in 0..nvtxs {
        if match_[i as usize] == UNMATCHED
            && ((xadj[(i + 1) as usize] - xadj[i as usize]) as usize) < maxdegree
        {
            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                colptr[adjncy[j as usize] as usize] += 1;
            }
        }
    }
    // MAKECSR(i, nvtxs, colptr)
    for i in 1..nvtxs {
        colptr[i as usize] += colptr[(i - 1) as usize];
    }
    for i in (1..=nvtxs).rev() {
        colptr[i as usize] = colptr[(i - 1) as usize];
    }
    colptr[0] = 0;

    let mut rowind = iwspacemalloc(ctrl, colptr[nvtxs as usize]);
    for pi in 0..nvtxs {
        let i = perm[pi as usize];
        if match_[i as usize] == UNMATCHED
            && ((xadj[(i + 1) as usize] - xadj[i as usize]) as usize) < maxdegree
        {
            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                let c = adjncy[j as usize] as usize;
                rowind[colptr[c] as usize] = i;
                colptr[c] += 1;
            }
        }
    }
    // SHIFTCSR(i, nvtxs, colptr)
    for i in (1..=nvtxs).rev() {
        colptr[i as usize] = colptr[(i - 1) as usize];
    }
    colptr[0] = 0;

    // compute matchings by going down the inverted index
    for pi in 0..nvtxs {
        let i = perm[pi as usize];
        if colptr[(i + 1) as usize] - colptr[i as usize] < 2 {
            continue;
        }

        let mut jj = colptr[(i + 1) as usize];
        let mut j = colptr[i as usize];
        while j < jj {
            if match_[rowind[j as usize] as usize] == UNMATCHED {
                jj -= 1;
                while jj > j {
                    if match_[rowind[jj as usize] as usize] == UNMATCHED {
                        let v = cnvtxs;
                        cnvtxs += 1;
                        cmap[rowind[jj as usize] as usize] = v;
                        cmap[rowind[j as usize] as usize] = v;
                        match_[rowind[j as usize] as usize] = rowind[jj as usize];
                        match_[rowind[jj as usize] as usize] = rowind[j as usize];
                        nunmatched = nunmatched.wrapping_sub(2);
                        break;
                    }
                    jj -= 1;
                }
            }
            j += 1;
        }
    }
    wspacepop(ctrl);

    *r_nunmatched = nunmatched;
    cnvtxs
}

/// HiGHS: coarsen.c::Match_2HopAll (lines 517-587)
pub fn Match_2HopAll(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    perm: &mut [idx_t],
    match_: &mut [idx_t],
    cnvtxs: idx_t,
    r_nunmatched: &mut usize,
    maxdegree: usize,
) -> idx_t {
    let mut cnvtxs = cnvtxs;
    let nvtxs = graph.nvtxs;
    let xadj: &[idx_t] = match &graph.xadj_orig {
        Some(x) => x.as_slice(),
        None => graph.xadj.as_slice(),
    };
    let adjncy: &[idx_t] = match &graph.adjncy_orig {
        Some(a) => a.as_slice(),
        None => graph.adjncy.as_slice(),
    };
    let cmap = &mut graph.cmap;

    let mut nunmatched = *r_nunmatched;
    let mask: idx_t = (IDX_MAX as usize / maxdegree) as idx_t;

    wspacepush(ctrl);

    // collapse vertices with identical adjacency lists
    let mut keys: Vec<ikv_t> = ikvwspacemalloc(ctrl, nunmatched as idx_t);
    let mut ncand: usize = 0;
    for pi in 0..nvtxs {
        let i = perm[pi as usize];
        let idegree = xadj[(i + 1) as usize] - xadj[i as usize];
        if match_[i as usize] == UNMATCHED && idegree > 1 && (idegree as usize) < maxdegree {
            let mut k: idx_t = 0;
            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                k += adjncy[j as usize] % mask;
            }
            keys[ncand].val = i;
            keys[ncand].key = ((((k % mask) as usize).wrapping_mul(maxdegree))
                .wrapping_add(idegree as usize)) as idx_t;
            ncand += 1;
        }
    }
    ikvsorti(ncand, &mut keys);

    let mut mark = iwspacemalloc(ctrl, nvtxs);
    iset(nvtxs as usize, 0, &mut mark);
    for pi in 0..ncand {
        let i = keys[pi].val;
        if match_[i as usize] != UNMATCHED {
            continue;
        }

        for j in xadj[i as usize]..xadj[(i + 1) as usize] {
            mark[adjncy[j as usize] as usize] = i;
        }

        for pk in (pi + 1)..ncand {
            let k = keys[pk].val;
            if match_[k as usize] != UNMATCHED {
                continue;
            }

            if keys[pi].key != keys[pk].key {
                break;
            }
            if xadj[(i + 1) as usize] - xadj[i as usize]
                != xadj[(k + 1) as usize] - xadj[k as usize]
            {
                break;
            }

            let mut jj = xadj[k as usize];
            while jj < xadj[(k + 1) as usize] {
                if mark[adjncy[jj as usize] as usize] != i {
                    break;
                }
                jj += 1;
            }
            if jj == xadj[(k + 1) as usize] {
                let v = cnvtxs;
                cnvtxs += 1;
                cmap[i as usize] = v;
                cmap[k as usize] = v;
                match_[i as usize] = k;
                match_[k as usize] = i;
                nunmatched = nunmatched.wrapping_sub(2);
                break;
            }
        }
    }
    wspacepop(ctrl);

    *r_nunmatched = nunmatched;
    cnvtxs
}

/// HiGHS: coarsen.c::PrintCGraphStats (lines 593-603)
/// Debug printer (`METIS_DBG_COARSEN`); never reached from HiPO.
pub fn PrintCGraphStats(ctrl: &ctrl_t, graph: &graph_t) {
    // HIGHS_ORDERING_PRINT expands to nothing under NDEBUG (oracle build).
    // C lines 593-603: all output is debug-only.
    let _ = ctrl;
    let _ = graph;
}

/// HiGHS: coarsen.c::CreateCoarseGraph (lines 612-846)
pub fn CreateCoarseGraph(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    cnvtxs: idx_t,
    match_: &[idx_t],
) {
    wspacepush(ctrl);

    let dovsize: i32 = 0;
    let dropedges = ctrl.dropedges;

    let mask: idx_t = HTLENGTH;

    let nvtxs = pool.get(graph).nvtxs;
    let ncon = pool.get(graph).ncon;

    let mut cnvtxs = cnvtxs;
    let mut keys: Vec<idx_t> = Vec::new();
    let mut noise: Vec<idx_t> = Vec::new();
    let mut medianewgts: Vec<idx_t> = Vec::new();

    // Setup structures for dropedges
    if dropedges != 0 {
        let mut nkeys: idx_t = 0;
        {
            let xadj = pool.get(graph).xadj_eff();
            for v in 0..nvtxs {
                nkeys = cpp_max(nkeys, xadj[(v + 1) as usize] - xadj[v as usize]);
            }
        }
        nkeys = 2 * nkeys + 1;

        keys = iwspacemalloc(ctrl, nkeys);
        noise = iwspacemalloc(ctrl, cnvtxs);
        medianewgts = iwspacemalloc(ctrl, cnvtxs);
        iset(cnvtxs as usize, -1, &mut medianewgts);

        for v in 0..cnvtxs {
            noise[v as usize] = irandInRange(128, &mut ctrl.rng_state);
        }
    }

    // Initialize the coarser graph
    let cgid = SetupCoarseGraph(pool, graph, cnvtxs, dovsize);
    let (g, cg) = pool.get2_mut(graph, cgid);
    let vwgt = &g.vwgt;
    let vsize = &g.vsize;
    let adjwgt = &g.adjwgt;
    let xadj = g.xadj_eff();
    let adjncy = g.adjncy_eff();
    let cmap = &g.cmap;

    let mut htable = iwspacemalloc(ctrl, mask + 1);
    iset((mask + 1) as usize, -1, &mut htable); // hash table
    let mut dtable = iwspacemalloc(ctrl, cnvtxs);
    iset(cnvtxs as usize, -1, &mut dtable); // direct table

    cg.xadj[0] = 0;
    cnvtxs = 0;
    let mut cnedges: idx_t = 0;
    for v in 0..nvtxs {
        let u = match_[v as usize];
        if u < v {
            continue;
        }

        // take care of the vertices
        if ncon == 1 {
            cg.vwgt[cnvtxs as usize] = vwgt[v as usize];
        } else {
            icopy(
                ncon as usize,
                &vwgt[(v * ncon) as usize..],
                &mut cg.vwgt[(cnvtxs * ncon) as usize..],
            );
        }

        if dovsize != 0 {
            cg.vsize[cnvtxs as usize] = vsize[v as usize];
        }

        if v != u {
            if ncon == 1 {
                cg.vwgt[cnvtxs as usize] += vwgt[u as usize];
            } else {
                iaxpy(
                    ncon as usize,
                    1,
                    &vwgt[(u * ncon) as usize..],
                    1,
                    &mut cg.vwgt[(cnvtxs * ncon) as usize..],
                    1,
                );
            }

            if dovsize != 0 {
                cg.vsize[cnvtxs as usize] += vsize[u as usize];
            }
        }

        // take care of the edges; the C pointers cadjncy / cadjwgt start at cnedges
        let base = cnedges as usize;
        let cadjncy = &mut cg.adjncy[base..];
        let cadjwgt = &mut cg.adjwgt[base..];
        let mut nedges: idx_t;

        if (xadj[(v + 1) as usize] - xadj[v as usize] + xadj[(u + 1) as usize] - xadj[u as usize])
            < (mask >> 2)
        {
            // use mask
            // put the ID of the contracted node itself at the start, so that it can be
            // removed easily
            htable[(cnvtxs & mask) as usize] = 0;
            cadjncy[0] = cnvtxs;
            nedges = 1;

            for j in xadj[v as usize]..xadj[(v + 1) as usize] {
                let k = cmap[adjncy[j as usize] as usize];
                let mut kk = k & mask;
                while htable[kk as usize] != -1 && cadjncy[htable[kk as usize] as usize] != k {
                    kk = (kk + 1) & mask;
                }
                let m = htable[kk as usize];
                if m == -1 {
                    cadjncy[nedges as usize] = k;
                    cadjwgt[nedges as usize] = adjwgt[j as usize];
                    htable[kk as usize] = nedges;
                    nedges += 1;
                } else {
                    cadjwgt[m as usize] += adjwgt[j as usize];
                }
            }

            if v != u {
                for j in xadj[u as usize]..xadj[(u + 1) as usize] {
                    let k = cmap[adjncy[j as usize] as usize];
                    let mut kk = k & mask;
                    while htable[kk as usize] != -1 && cadjncy[htable[kk as usize] as usize] != k {
                        kk = (kk + 1) & mask;
                    }
                    let m = htable[kk as usize];
                    if m == -1 {
                        cadjncy[nedges as usize] = k;
                        cadjwgt[nedges as usize] = adjwgt[j as usize];
                        htable[kk as usize] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m as usize] += adjwgt[j as usize];
                    }
                }
            }

            // reset the htable -- reverse order (LIFO) is critical to prevent cadjncy[-1]
            // indexing due to a remove of an earlier entry
            let mut j = nedges - 1;
            while j >= 0 {
                let k = cadjncy[j as usize];
                let mut kk = k & mask;
                while cadjncy[htable[kk as usize] as usize] != k {
                    kk = (kk + 1) & mask;
                }
                htable[kk as usize] = -1;
                j -= 1;
            }

            // remove the contracted vertex from the list
            nedges -= 1;
            cadjncy[0] = cadjncy[nedges as usize];
            cadjwgt[0] = cadjwgt[nedges as usize];
        } else {
            nedges = 0;
            for j in xadj[v as usize]..xadj[(v + 1) as usize] {
                let k = cmap[adjncy[j as usize] as usize];
                let m = dtable[k as usize];
                if m == -1 {
                    cadjncy[nedges as usize] = k;
                    cadjwgt[nedges as usize] = adjwgt[j as usize];
                    dtable[k as usize] = nedges;
                    nedges += 1;
                } else {
                    cadjwgt[m as usize] += adjwgt[j as usize];
                }
            }

            if v != u {
                for j in xadj[u as usize]..xadj[(u + 1) as usize] {
                    let k = cmap[adjncy[j as usize] as usize];
                    let m = dtable[k as usize];
                    if m == -1 {
                        cadjncy[nedges as usize] = k;
                        cadjwgt[nedges as usize] = adjwgt[j as usize];
                        dtable[k as usize] = nedges;
                        nedges += 1;
                    } else {
                        cadjwgt[m as usize] += adjwgt[j as usize];
                    }
                }

                // Remove the contracted self-loop, when present
                let j = dtable[cnvtxs as usize];
                if j != -1 {
                    nedges -= 1;
                    cadjncy[j as usize] = cadjncy[nedges as usize];
                    cadjwgt[j as usize] = cadjwgt[nedges as usize];
                    dtable[cnvtxs as usize] = -1;
                }
            }

            // Zero out the dtable
            for j in 0..nedges {
                dtable[cadjncy[j as usize] as usize] = -1;
            }
        }

        // Determine the median weight of the incident edges, which will be used
        // to keep an edge (u, v) iff wgt(u, v) >= min(medianewgts[u], medianewgts[v])
        if dropedges != 0 {
            medianewgts[cnvtxs as usize] = 8; // default for island nodes
            if nedges > 0 {
                for j in 0..nedges {
                    keys[j as usize] = (cadjwgt[j as usize] << 8)
                        + noise[cnvtxs as usize]
                        + noise[cadjncy[j as usize] as usize];
                }
                isortd(nedges as usize, &mut keys);
                medianewgts[cnvtxs as usize] = keys[cpp_min(
                    nedges - 1,
                    (xadj[(v + 1) as usize] - xadj[v as usize] + xadj[(u + 1) as usize]
                        - xadj[u as usize])
                        >> 1,
                ) as usize];
            }
        }

        cnedges += nedges;
        cnvtxs += 1;
        cg.xadj[cnvtxs as usize] = cnedges;
    }

    // compact the adjacency structure of the coarser graph to keep only +ve edges
    if dropedges != 0 {
        let mut droppedewgt: idx_t = 0;

        cnedges = 0;
        for u in 0..cnvtxs {
            let istart = cg.xadj[u as usize];
            let iend = cg.xadj[(u + 1) as usize];
            for j in istart..iend {
                let v = cg.adjncy[j as usize];
                if (cg.adjwgt[j as usize] << 8) + noise[u as usize] + noise[v as usize]
                    >= cpp_min(medianewgts[u as usize], medianewgts[v as usize])
                {
                    cg.adjncy[cnedges as usize] = cg.adjncy[j as usize];
                    cg.adjwgt[cnedges as usize] = cg.adjwgt[j as usize];
                    cnedges += 1;
                } else {
                    droppedewgt += cg.adjwgt[j as usize];
                }
            }
            cg.xadj[u as usize] = cnedges;
        }
        // SHIFTCSR(j, cnvtxs, cxadj)
        for j in (1..=cnvtxs).rev() {
            cg.xadj[j as usize] = cg.xadj[(j - 1) as usize];
        }
        cg.xadj[0] = 0;

        cg.droppedewgt = droppedewgt;
    }

    cg.nedges = cnedges;

    for j in 0..ncon {
        let tv = isum(cg.nvtxs as usize, &cg.vwgt[j as usize..], ncon as usize);
        cg.tvwgt[j as usize] = tv;
        let denom = if tv > 0 { tv } else { 1 };
        cg.invtvwgt[j as usize] = (1.0f64 / (denom as f64)) as real_t;
    }

    ReAdjustMemory(ctrl, g, cg);

    wspacepop(ctrl);
}

/// HiGHS: coarsen.c::SetupCoarseGraph (lines 853-881)
/// Allocates the coarse graph in the pool and links `coarser` / `finer`.
pub fn SetupCoarseGraph(
    pool: &mut GraphPool,
    graph: GraphId,
    cnvtxs: idx_t,
    dovsize: i32,
) -> GraphId {
    let mut cgraph = CreateGraph();

    cgraph.nvtxs = cnvtxs;
    cgraph.ncon = pool.get(graph).ncon;

    cgraph.finer = Some(graph);

    // NOTE: The +1 in the adjwgt/adjncy is to allow the optimization of self-loop
    // detection by adding ahead of time the self-loop.
    let nedges = pool.get(graph).nedges;
    cgraph.xadj = imalloc((cnvtxs + 1) as usize);
    cgraph.adjncy = imalloc((nedges + 1) as usize);
    cgraph.adjwgt = imalloc((nedges + 1) as usize);
    cgraph.vwgt = imalloc((cgraph.ncon * cnvtxs) as usize);
    cgraph.tvwgt = imalloc(cgraph.ncon as usize);
    cgraph.invtvwgt = rmalloc(cgraph.ncon as usize);

    if dovsize != 0 {
        cgraph.vsize = imalloc(cnvtxs as usize);
    }

    let cid = pool.alloc(cgraph);
    pool.get_mut(graph).coarser = Some(cid);
    cid
}

/// HiGHS: coarsen.c::ReAdjustMemory (lines 889-895)
pub fn ReAdjustMemory(ctrl: &mut ctrl_t, graph: &graph_t, cgraph: &mut graph_t) {
    let _ = ctrl;
    if cgraph.nedges > 10000 && (cgraph.nedges as f64) < 0.9 * (graph.nedges as f64) {
        cgraph.adjncy.resize(cgraph.nedges as usize, 0);
        cgraph.adjwgt.resize(cgraph.nedges as usize, 0);
    }
}

#[cfg(test)]
mod tests {
    //! Hand-derived sanity checks, traced by hand through
    //! coarsen.c; not oracle fixtures. Paths that reach the LGPL-held sorts
    //! (`ikvsorti`, `isortd`) are not exercised.
    use super::*;

    fn graph_from(xadj: Vec<idx_t>, adjncy: Vec<idx_t>) -> graph_t {
        let mut g = CreateGraph();
        g.nvtxs = (xadj.len() - 1) as idx_t;
        g.nedges = adjncy.len() as idx_t;
        g.ncon = 1;
        g.vwgt = vec![1; g.nvtxs as usize];
        g.adjwgt = vec![1; adjncy.len()];
        g.cmap = vec![0; g.nvtxs as usize];
        g.xadj = xadj;
        g.adjncy = adjncy;
        g
    }

    #[test]
    fn match_2hop_any_pairs_leaves_of_path() {
        // Path 0-1-2: 0 and 2 share neighbour 1 (coarsen.c 444-506).
        let mut ctrl = ctrl_t::default();
        let mut g = graph_from(vec![0, 1, 3, 4], vec![1, 0, 2, 1]);
        let mut perm = vec![0, 1, 2];
        let mut match_ = vec![UNMATCHED; 3];
        let mut nunmatched = 3usize;
        let cnvtxs = Match_2HopAny(
            &mut ctrl,
            &mut g,
            &mut perm,
            &mut match_,
            0,
            &mut nunmatched,
            2,
        );
        assert_eq!(cnvtxs, 1);
        assert_eq!(nunmatched, 1);
        assert_eq!(match_, vec![2, UNMATCHED, 0]);
        assert_eq!(g.cmap[0], 0);
        assert_eq!(g.cmap[2], 0);
        assert_eq!(ctrl.mcore, 0);
    }

    #[test]
    fn create_coarse_graph_hash_path() {
        // Path 0-1-2-3 contracted as {0,1}, {2,3} (coarsen.c 612-846, mask path).
        let mut ctrl = ctrl_t::default();
        let mut g = graph_from(vec![0, 1, 3, 5, 6], vec![1, 0, 2, 1, 3, 2]);
        g.cmap = vec![0, 0, 1, 1];
        let mut pool = GraphPool::new();
        let gid = pool.alloc(g);
        CreateCoarseGraph(&mut ctrl, &mut pool, gid, 2, &[1, 0, 3, 2]);
        let cid = pool.get(gid).coarser.expect("coarser linked");
        let cg = pool.get(cid);
        assert_eq!(cg.finer, Some(gid));
        assert_eq!(cg.nvtxs, 2);
        assert_eq!(cg.nedges, 2);
        assert_eq!(&cg.xadj[..], &[0, 1, 2]);
        assert_eq!(&cg.adjncy[..2], &[1, 0]);
        assert_eq!(&cg.adjwgt[..2], &[1, 1]);
        assert_eq!(cg.vwgt, vec![2, 2]);
        assert_eq!(cg.tvwgt, vec![4]);
        assert_eq!(cg.invtvwgt, vec![0.25]);
        assert_eq!(ctrl.mcore, 0);
    }

    #[test]
    fn create_coarse_graph_direct_path() {
        // Star with 1024 leaves, nothing matched: degree sum 2048 >= HTLENGTH>>2
        // sends the centre through the dtable path.
        let nleaves: idx_t = 1024;
        let mut xadj = vec![0, nleaves];
        let mut adjncy: Vec<idx_t> = (1..=nleaves).collect();
        for _ in 0..nleaves {
            adjncy.push(0);
            xadj.push(adjncy.len() as idx_t);
        }
        let mut ctrl = ctrl_t::default();
        let mut g = graph_from(xadj, adjncy);
        let n = g.nvtxs;
        g.cmap = (0..n).collect();
        let match_: Vec<idx_t> = (0..n).collect();
        let mut pool = GraphPool::new();
        let gid = pool.alloc(g);
        CreateCoarseGraph(&mut ctrl, &mut pool, gid, n, &match_);
        let cid = pool.get(gid).coarser.expect("coarser linked");
        let cg = pool.get(cid);
        assert_eq!(cg.nedges, 2 * nleaves);
        assert_eq!(cg.xadj[1], nleaves);
        let expect: Vec<idx_t> = (1..=nleaves).collect();
        assert_eq!(&cg.adjncy[..nleaves as usize], &expect[..]);
        assert_eq!(cg.adjncy[nleaves as usize], 0);
        assert_eq!(cg.tvwgt, vec![n]);
    }
}
