//! HiGHS source: extern/metis/libmetis/initpart.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Initial separators of the coarsest graph: `InitSeparator` dispatches on
//! `ctrl->iptype` (edge: random / grow bisection then `ConstructSeparator`;
//! node: `GrowBisectionNode`). Each trial draws from `ctrl->rng_state`.

use super::balance::Balance2Way;
use super::compress::ConstructSeparator;
use super::fm::FM_2WayRefine;
use super::gklib::{gk_errexit, icopy, imalloc, irandArrayPermute, irandInRange, iset};
use super::options::Setup2WayBalMultipliers;
use super::refine::{
    Allocate2WayPartitionMemory, Compute2WayNodePartitionParams, Compute2WayPartitionParams,
};
use super::sfm::{FM_2WayNodeRefine1Sided, FM_2WayNodeRefine2Sided};
use super::wspace::{iwspacemalloc, wspacepop, wspacepush};
use super::{ctrl_t, graph_t, idx_t, real_t};

/// HiGHS: initpart.c::InitSeparator (lines 19-54)
pub fn InitSeparator(ctrl: &mut ctrl_t, graph: &mut graph_t, niparts: idx_t) {
    let mut ntpwgts: [real_t; 2] = [0.5, 0.5];
    let dbglvl = ctrl.dbglvl;

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) dropped (oracle config dbglvl == 0)
    // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...) dropped (oracle config dbglvl == 0)

    // this is required for the cut-based part of the refinement
    Setup2WayBalMultipliers(ctrl, graph, &ntpwgts);

    match ctrl.iptype {
        super::METIS_IPTYPE_EDGE => {
            if graph.nedges == 0 {
                RandomBisection(ctrl, graph, &ntpwgts, niparts);
            } else {
                GrowBisection(ctrl, graph, &ntpwgts, niparts);
            }

            Compute2WayPartitionParams(ctrl, graph);
            ConstructSeparator(ctrl, graph);
        }

        super::METIS_IPTYPE_NODE => {
            GrowBisectionNode(ctrl, graph, &ntpwgts, niparts);
        }

        _ => {
            gk_errexit(&format!("Unknown iptype of {}\n", ctrl.iptype));
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_IPART, ...) dropped (oracle config dbglvl == 0)

    ctrl.dbglvl = dbglvl;
}

/// HiGHS: initpart.c::RandomBisection (lines 63-127)
pub fn RandomBisection(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t], niparts: idx_t) {
    let mut i: idx_t;
    let mut ii: idx_t;
    let mut j: idx_t;
    let mut k: idx_t;
    let mut pwgts: [idx_t; 2] = [0; 2];
    let mut from: idx_t;
    let mut me: idx_t;
    let mut bestcut: idx_t = 0;
    let mut icut: idx_t;
    let mut mincut: idx_t;
    let mut inbfs: idx_t;

    wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    let xadj = graph.xadj_eff().to_vec();
    let vwgt = graph.vwgt.clone();
    let adjncy = graph.adjncy_eff().to_vec();
    let adjwgt = graph.adjwgt.clone();

    Allocate2WayPartitionMemory(ctrl, graph);

    let mut bestwhere = iwspacemalloc(ctrl, nvtxs);
    let mut perm = iwspacemalloc(ctrl, nvtxs);

    // zeromaxpwgt = ctrl->ubfactors[0]*graph->tvwgt[0]*ntpwgts[0];
    // f32 * (idx_t->f32) * f32, then truncated to idx_t
    let zeromaxpwgt: idx_t =
        ((ctrl.ubfactors[0] * (graph.tvwgt[0] as real_t)) * ntpwgts[0]) as idx_t;

    inbfs = 0;
    while inbfs < niparts {
        iset(nvtxs as usize, 1, &mut graph.where_);

        if inbfs > 0 {
            irandArrayPermute(nvtxs, &mut perm, nvtxs / 2, 1, &mut ctrl.rng_state);
            pwgts[1] = graph.tvwgt[0];
            pwgts[0] = 0;

            ii = 0;
            while ii < nvtxs {
                i = perm[ii as usize];
                if pwgts[0] + vwgt[i as usize] < zeromaxpwgt {
                    graph.where_[i as usize] = 0;
                    pwgts[0] += vwgt[i as usize];
                    pwgts[1] -= vwgt[i as usize];
                    if pwgts[0] > zeromaxpwgt {
                        break;
                    }
                }
                ii += 1;
            }
        }

        // Do some partition refinement
        Compute2WayPartitionParams(ctrl, graph);

        Balance2Way(ctrl, graph, ntpwgts);

        FM_2WayRefine(ctrl, graph, ntpwgts, 4);

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            icopy(nvtxs as usize, &graph.where_, &mut bestwhere);
            if bestcut == 0 {
                break;
            }
        }

        inbfs += 1;
    }

    graph.mincut = bestcut;
    icopy(nvtxs as usize, &bestwhere, &mut graph.where_);

    wspacepop(ctrl);
}

/// HiGHS: initpart.c::GrowBisection (lines 136-252)
pub fn GrowBisection(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t], niparts: idx_t) {
    let mut i: idx_t;
    let mut j: idx_t;
    let mut k: idx_t;
    let mut drain: idx_t;
    let mut nleft: idx_t;
    let mut first: idx_t;
    let mut last: idx_t;
    let mut pwgts: [idx_t; 2] = [0; 2];
    let mut from: idx_t;
    let mut me: idx_t;
    let mut bestcut: idx_t = 0;
    let mut icut: idx_t;
    let mut mincut: idx_t;
    let mut inbfs: idx_t;

    wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    let xadj = graph.xadj_eff().to_vec();
    let vwgt = graph.vwgt.clone();
    let adjncy = graph.adjncy_eff().to_vec();
    let adjwgt = graph.adjwgt.clone();

    Allocate2WayPartitionMemory(ctrl, graph);

    let mut bestwhere = iwspacemalloc(ctrl, nvtxs);
    let mut queue = iwspacemalloc(ctrl, nvtxs);
    let mut touched = iwspacemalloc(ctrl, nvtxs);

    // onemaxpwgt = ctrl->ubfactors[0]*graph->tvwgt[0]*ntpwgts[1];
    // f32 * (idx_t->f32) * f32, then truncated to idx_t
    let onemaxpwgt: idx_t =
        ((ctrl.ubfactors[0] * (graph.tvwgt[0] as real_t)) * ntpwgts[1]) as idx_t;

    // oneminpwgt = (1.0/ctrl->ubfactors[0])*graph->tvwgt[0]*ntpwgts[1];
    // 1.0 is double, so computed in f64 then truncated to idx_t
    let oneminpwgt: idx_t = ((1.0 / (ctrl.ubfactors[0] as f64))
        * (graph.tvwgt[0] as f64)
        * (ntpwgts[1] as f64)) as idx_t;

    inbfs = 0;
    while inbfs < niparts {
        iset(nvtxs as usize, 1, &mut graph.where_);

        iset(nvtxs as usize, 0, &mut touched);

        pwgts[1] = graph.tvwgt[0];
        pwgts[0] = 0;

        queue[0] = irandInRange(nvtxs, &mut ctrl.rng_state);
        touched[queue[0] as usize] = 1;
        first = 0;
        last = 1;
        nleft = nvtxs - 1;
        drain = 0;

        // Start the BFS from queue to get a partition
        loop {
            if first == last {
                // Empty. Disconnected graph!
                if nleft == 0 || drain != 0 {
                    break;
                }

                k = irandInRange(nleft, &mut ctrl.rng_state);
                i = 0;
                while i < nvtxs {
                    if touched[i as usize] == 0 {
                        if k == 0 {
                            break;
                        } else {
                            k -= 1;
                        }
                    }
                    i += 1;
                }

                queue[0] = i;
                touched[i as usize] = 1;
                first = 0;
                last = 1;
                nleft -= 1;
            }

            i = queue[first as usize];
            first += 1;

            if pwgts[0] > 0 && pwgts[1] - vwgt[i as usize] < oneminpwgt {
                drain = 1;
                continue;
            }

            graph.where_[i as usize] = 0;
            // INC_DEC(pwgts[0], pwgts[1], vwgt[i])
            pwgts[0] += vwgt[i as usize];
            pwgts[1] -= vwgt[i as usize];

            if pwgts[1] <= onemaxpwgt {
                break;
            }

            drain = 0;
            j = xadj[i as usize];
            while j < xadj[(i + 1) as usize] {
                k = adjncy[j as usize];
                if touched[k as usize] == 0 {
                    queue[last as usize] = k;
                    last += 1;
                    touched[k as usize] = 1;
                    nleft -= 1;
                }
                j += 1;
            }
        }

        // Check to see if we hit any bad limiting cases
        if pwgts[1] == 0 {
            let r = irandInRange(nvtxs, &mut ctrl.rng_state) as usize;
            graph.where_[r] = 1;
        }
        if pwgts[0] == 0 {
            let r = irandInRange(nvtxs, &mut ctrl.rng_state) as usize;
            graph.where_[r] = 0;
        }

        // Do some partition refinement
        Compute2WayPartitionParams(ctrl, graph);

        Balance2Way(ctrl, graph, ntpwgts);

        FM_2WayRefine(ctrl, graph, ntpwgts, ctrl.niter);

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            icopy(nvtxs as usize, &graph.where_, &mut bestwhere);
            if bestcut == 0 {
                break;
            }
        }

        inbfs += 1;
    }

    graph.mincut = bestcut;
    icopy(nvtxs as usize, &bestwhere, &mut graph.where_);

    wspacepop(ctrl);
}

/// HiGHS: initpart.c::GrowBisectionNode (lines 261-385)
pub fn GrowBisectionNode(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    ntpwgts: &[real_t],
    niparts: idx_t,
) {
    let mut i: idx_t;
    let mut j: idx_t;
    let mut k: idx_t;
    let mut drain: idx_t;
    let mut nleft: idx_t;
    let mut first: idx_t;
    let mut last: idx_t;
    let mut pwgts: [idx_t; 2] = [0; 2];
    let mut from: idx_t;
    let mut me: idx_t;
    let mut bestcut: idx_t = 0;
    let mut icut: idx_t;
    let mut mincut: idx_t;
    let mut inbfs: idx_t;

    wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    let xadj = graph.xadj_eff().to_vec();
    let vwgt = graph.vwgt.clone();
    let adjncy = graph.adjncy_eff().to_vec();
    let adjwgt = graph.adjwgt.clone();

    let mut bestwhere = iwspacemalloc(ctrl, nvtxs);
    let mut queue = iwspacemalloc(ctrl, nvtxs);
    let mut touched = iwspacemalloc(ctrl, nvtxs);

    // onemaxpwgt = ctrl->ubfactors[0]*graph->tvwgt[0]*0.5;
    // f32 * (idx_t->f32) in f32; 0.5 is a double literal, so the f32 product is
    // promoted and the last multiply is f64, then truncated to idx_t
    let onemaxpwgt: idx_t =
        (((ctrl.ubfactors[0] * (graph.tvwgt[0] as real_t)) as f64) * 0.5f64) as idx_t;

    // oneminpwgt = (1.0/ctrl->ubfactors[0])*graph->tvwgt[0]*0.5;
    // 1.0 is double, so computed in f64 then truncated to idx_t
    let oneminpwgt: idx_t =
        ((1.0 / (ctrl.ubfactors[0] as f64)) * (graph.tvwgt[0] as f64) * 0.5f64) as idx_t;

    // Allocate refinement memory. Allocate sufficient memory for both edge and node
    graph.pwgts = imalloc(3);
    graph.where_ = imalloc(nvtxs as usize);
    graph.bndptr = imalloc(nvtxs as usize);
    graph.bndind = imalloc(nvtxs as usize);
    graph.id = imalloc(nvtxs as usize);
    graph.ed = imalloc(nvtxs as usize);
    graph.nrinfo = vec![super::nrinfo_t::default(); nvtxs as usize];

    inbfs = 0;
    while inbfs < niparts {
        iset(nvtxs as usize, 1, &mut graph.where_);
        iset(nvtxs as usize, 0, &mut touched);

        pwgts[1] = graph.tvwgt[0];
        pwgts[0] = 0;

        queue[0] = irandInRange(nvtxs, &mut ctrl.rng_state);
        touched[queue[0] as usize] = 1;
        first = 0;
        last = 1;
        nleft = nvtxs - 1;
        drain = 0;

        // Start the BFS from queue to get a partition
        loop {
            if first == last {
                // Empty. Disconnected graph!
                if nleft == 0 || drain != 0 {
                    break;
                }

                k = irandInRange(nleft, &mut ctrl.rng_state);
                i = 0;
                while i < nvtxs {
                    // select the kth untouched vertex
                    if touched[i as usize] == 0 {
                        if k == 0 {
                            break;
                        } else {
                            k -= 1;
                        }
                    }
                    i += 1;
                }

                queue[0] = i;
                touched[i as usize] = 1;
                first = 0;
                last = 1;
                nleft -= 1;
            }

            i = queue[first as usize];
            first += 1;

            if pwgts[1] - vwgt[i as usize] < oneminpwgt {
                drain = 1;
                continue;
            }

            graph.where_[i as usize] = 0;
            // INC_DEC(pwgts[0], pwgts[1], vwgt[i])
            pwgts[0] += vwgt[i as usize];
            pwgts[1] -= vwgt[i as usize];

            if pwgts[1] <= onemaxpwgt {
                break;
            }

            drain = 0;
            j = xadj[i as usize];
            while j < xadj[(i + 1) as usize] {
                k = adjncy[j as usize];
                if touched[k as usize] == 0 {
                    queue[last as usize] = k;
                    last += 1;
                    touched[k as usize] = 1;
                    nleft -= 1;
                }
                j += 1;
            }
        }

        // Do some partition refinement
        Compute2WayPartitionParams(ctrl, graph);
        Balance2Way(ctrl, graph, ntpwgts);
        FM_2WayRefine(ctrl, graph, ntpwgts, 4);

        // Construct and refine the vertex separator
        i = 0;
        while i < graph.nbnd {
            j = graph.bndind[i as usize];
            if xadj[(j + 1) as usize] - xadj[j as usize] > 0 {
                // ignore islands
                graph.where_[j as usize] = 2;
            }
            i += 1;
        }

        Compute2WayNodePartitionParams(ctrl, graph);
        FM_2WayNodeRefine2Sided(ctrl, graph, 1);
        FM_2WayNodeRefine1Sided(ctrl, graph, 4);

        if inbfs == 0 || bestcut > graph.mincut {
            bestcut = graph.mincut;
            icopy(nvtxs as usize, &graph.where_, &mut bestwhere);
        }

        inbfs += 1;
    }

    graph.mincut = bestcut;
    icopy(nvtxs as usize, &bestwhere, &mut graph.where_);

    wspacepop(ctrl);
}

#[cfg(test)]
mod tests {
    //! Unit tests. Hand-derived invariants of initpart.c:
    //! edge bisections leave `where` in {0,1} with `mincut` the cut of the
    //! returned partition (lines 117-126 / 235-251); node bisections leave a
    //! vertex separator (`where == 2`) with `mincut` its weight (lines 362-384).
    use super::*;
    use crate::ipm::hipo::orderings::metis::graph::SetupGraph;
    use crate::ipm::hipo::orderings::metis::options::SetupCtrl;
    use crate::ipm::hipo::orderings::metis::{METIS_IPTYPE_EDGE, METIS_IPTYPE_NODE};
    use crate::ipm::hipo::orderings::metis::{
        METIS_NOPTIONS, METIS_OP_OMETIS, METIS_OPTION_IPTYPE,
    };

    fn setup(iptype: idx_t, xadj: &[idx_t], adjncy: &[idx_t]) -> (Box<ctrl_t>, graph_t) {
        let mut options = [-1 as idx_t; METIS_NOPTIONS];
        options[METIS_OPTION_IPTYPE] = iptype;
        // ometis.c: SetupCtrl(METIS_OP_OMETIS, options, 1, 3, NULL, NULL)
        let mut ctrl = SetupCtrl(METIS_OP_OMETIS, Some(&options), 1, 3, None, None).unwrap();
        let n = (xadj.len() - 1) as idx_t;
        let graph = SetupGraph(&mut ctrl, n, 1, xadj, adjncy, None, None, None);
        (ctrl, graph)
    }

    /// 4x4 grid, vertex r*4+c.
    fn grid4() -> (Vec<idx_t>, Vec<idx_t>) {
        let mut xadj = vec![0];
        let mut adjncy = Vec::new();
        for r in 0..4 {
            for c in 0..4 {
                if r > 0 {
                    adjncy.push((r - 1) * 4 + c);
                }
                if c > 0 {
                    adjncy.push(r * 4 + c - 1);
                }
                if c < 3 {
                    adjncy.push(r * 4 + c + 1);
                }
                if r < 3 {
                    adjncy.push((r + 1) * 4 + c);
                }
                xadj.push(adjncy.len() as idx_t);
            }
        }
        (xadj, adjncy)
    }

    /// Two disjoint triangles {0,1,2} and {3,4,5}.
    fn two_triangles() -> (Vec<idx_t>, Vec<idx_t>) {
        (
            vec![0, 2, 4, 6, 8, 10, 12],
            vec![1, 2, 0, 2, 0, 1, 4, 5, 3, 5, 3, 4],
        )
    }

    fn edge_cut(g: &graph_t) -> idx_t {
        let xadj = g.xadj_eff();
        let adjncy = g.adjncy_eff();
        let mut cut = 0;
        for i in 0..g.nvtxs as usize {
            for j in xadj[i] as usize..xadj[i + 1] as usize {
                if g.where_[i] != g.where_[adjncy[j] as usize] {
                    cut += g.adjwgt[j];
                }
            }
        }
        cut / 2
    }

    fn assert_separator(g: &graph_t) {
        let xadj = g.xadj_eff();
        let adjncy = g.adjncy_eff();
        let mut sepwgt = 0;
        for i in 0..g.nvtxs as usize {
            let wi = g.where_[i];
            assert!((0..=2).contains(&wi), "where[{i}] = {wi}");
            if wi == 2 {
                sepwgt += g.vwgt[i];
                continue;
            }
            for j in xadj[i] as usize..xadj[i + 1] as usize {
                let wk = g.where_[adjncy[j] as usize];
                assert!(wk == wi || wk == 2, "edge {i}-{} crosses 0/1", adjncy[j]);
            }
        }
        assert_eq!(g.mincut, sepwgt);
    }

    #[test]
    fn grow_bisection_grid_cut_matches() {
        let (xadj, adjncy) = grid4();
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
        let ntpwgts: [real_t; 2] = [0.5, 0.5];
        Setup2WayBalMultipliers(&mut ctrl, &g, &ntpwgts);
        GrowBisection(&mut ctrl, &mut g, &ntpwgts, 4);
        assert_eq!(ctrl.mcore, 0);
        assert!(g.where_.iter().all(|&w| w == 0 || w == 1));
        assert!(g.where_.contains(&0) && g.where_.contains(&1));
        assert_eq!(g.mincut, edge_cut(&g));
        // 4x4 grid: a straight cut of 4 edges exists and FM finds it.
        assert!(g.mincut <= 4, "mincut {}", g.mincut);
    }

    #[test]
    fn grow_bisection_disconnected_zero_cut() {
        // BFS empties the queue on the first triangle (lines 172-189), then
        // pwgts[1] - vwgt < oneminpwgt drains; a zero cut stops the trials.
        let (xadj, adjncy) = two_triangles();
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
        let ntpwgts: [real_t; 2] = [0.5, 0.5];
        Setup2WayBalMultipliers(&mut ctrl, &g, &ntpwgts);
        GrowBisection(&mut ctrl, &mut g, &ntpwgts, 7);
        assert_eq!(g.mincut, edge_cut(&g));
        assert_eq!(g.mincut, 0);
        let w0 = g.where_[0];
        assert!(g.where_[..3].iter().all(|&w| w == w0));
        assert!(g.where_[3..].iter().all(|&w| w == 1 - w0));
    }

    #[test]
    fn random_bisection_grid_cut_matches() {
        let (xadj, adjncy) = grid4();
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
        let ntpwgts: [real_t; 2] = [0.5, 0.5];
        Setup2WayBalMultipliers(&mut ctrl, &g, &ntpwgts);
        RandomBisection(&mut ctrl, &mut g, &ntpwgts, 4);
        assert_eq!(ctrl.mcore, 0);
        assert!(g.where_.iter().all(|&w| w == 0 || w == 1));
        assert_eq!(g.mincut, edge_cut(&g));
    }

    #[test]
    fn grow_bisection_node_grid_separator() {
        let (xadj, adjncy) = grid4();
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_NODE, &xadj, &adjncy);
        let ntpwgts: [real_t; 2] = [0.5, 0.5];
        Setup2WayBalMultipliers(&mut ctrl, &g, &ntpwgts);
        GrowBisectionNode(&mut ctrl, &mut g, &ntpwgts, 4);
        assert_eq!(ctrl.mcore, 0);
        assert_eq!(g.pwgts.len(), 3);
        assert_eq!(g.nrinfo.len(), 16);
        assert_separator(&g);
        assert!(g.mincut > 0);
    }

    #[test]
    fn init_separator_edge_and_node_grid() {
        for iptype in [METIS_IPTYPE_EDGE, METIS_IPTYPE_NODE] {
            let (xadj, adjncy) = grid4();
            let (mut ctrl, mut g) = setup(iptype, &xadj, &adjncy);
            let dbglvl = ctrl.dbglvl;
            InitSeparator(&mut ctrl, &mut g, 5);
            assert_eq!(ctrl.dbglvl, dbglvl);
            assert_eq!(ctrl.mcore, 0);
            assert_separator(&g);
        }
    }

    #[test]
    fn init_separator_edgeless_uses_random_bisection() {
        // nedges == 0 -> RandomBisection (line 33), then ConstructSeparator.
        let xadj = [0, 0, 0, 0, 0];
        let adjncy: [idx_t; 0] = [];
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
        assert_eq!(g.nedges, 0);
        InitSeparator(&mut ctrl, &mut g, 5);
        assert_separator(&g);
        assert_eq!(g.mincut, 0);
    }

    #[test]
    fn init_separator_deterministic_rng() {
        let (xadj, adjncy) = grid4();
        let run = || {
            let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
            InitSeparator(&mut ctrl, &mut g, 5);
            (g.where_.clone(), g.mincut, ctrl.rng_state)
        };
        assert_eq!(run(), run());
    }

    #[test]
    #[should_panic]
    fn init_separator_unknown_iptype_aborts() {
        let (xadj, adjncy) = grid4();
        let (mut ctrl, mut g) = setup(METIS_IPTYPE_EDGE, &xadj, &adjncy);
        ctrl.iptype = 7;
        InitSeparator(&mut ctrl, &mut g, 1);
    }
}
