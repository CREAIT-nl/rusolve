//! HiGHS source: extern/metis/libmetis/ometis.c (lines 307-661)
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Node bisection at one level (`L2` tries `nseps`-style multiple coarsenings,
//! `L1` one), graph splitting along the separator, and the MMD leaf ordering.

use super::super::coarsen::{CoarsenGraph, CoarsenGraphNlevels};
use super::super::graph::{FreeRData, SetupGraph_tvwgt, SetupSplitGraph};
use super::super::initpart::InitSeparator;
use super::super::mmd::genmmd;
use super::super::refine::Refine2WayNode;
use super::super::wspace;
use super::super::{GraphId, GraphPool, ctrl_t, graph_t, idx_t};

/// HiGHS: ometis.c::MlevelNodeBisectionL2 (lines 313-357)
pub fn MlevelNodeBisectionL2(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    niparts: idx_t,
) {
    // HiGHS: ometis.c line 315: idx_t i, mincut, nruns=5;
    let nruns: idx_t = 5;

    // HiGHS: ometis.c lines 319-323: if the graph is small, just find a single vertex separator
    if pool.get(graph).nvtxs < 5000 {
        MlevelNodeBisectionL1(ctrl, pool, graph, niparts);
        return;
    }

    // HiGHS: ometis.c line 325: WCOREPUSH;
    wspace::wspacepush(ctrl);

    // HiGHS: ometis.c line 327: ctrl->CoarsenTo = gk_max(100, graph->nvtxs/30);
    // gk_max is a C ternary: if a > b { a } else { b }
    {
        let a: idx_t = 100;
        let b: idx_t = pool.get(graph).nvtxs / 30;
        ctrl.CoarsenTo = if a > b { a } else { b };
    }

    // HiGHS: ometis.c line 329: cgraph = CoarsenGraphNlevels(ctrl, graph, 4);
    let cgraph = CoarsenGraphNlevels(ctrl, pool, graph, 4);

    // HiGHS: ometis.c line 331: bestwhere = iwspacemalloc(ctrl, cgraph->nvtxs);
    let cnvtxs = pool.get(cgraph).nvtxs;
    let mut bestwhere = wspace::iwspacemalloc(ctrl, cnvtxs);

    // HiGHS: ometis.c line 333: mincut = graph->tvwgt[0];
    let mut mincut = pool.get(graph).tvwgt[0];

    // HiGHS: ometis.c lines 334-348
    let mut i: idx_t = 0;
    while i < nruns {
        // HiGHS: ometis.c line 335: MlevelNodeBisectionL1(ctrl, cgraph, 0.7*niparts);
        // 0.7 is a double literal; 0.7*niparts promotes to double then truncates to idx_t.
        let niparts_arg: idx_t = (0.7 * (niparts as f64)) as idx_t;
        MlevelNodeBisectionL1(ctrl, pool, cgraph, niparts_arg);

        // HiGHS: ometis.c lines 337-341
        let cgraph_mincut = pool.get(cgraph).mincut;
        if i == 0 || cgraph_mincut < mincut {
            mincut = cgraph_mincut;
            if i < nruns - 1 {
                // HiGHS: ometis.c line 340: icopy(cgraph->nvtxs, cgraph->where, bestwhere);
                let cg = pool.get(cgraph);
                for k in 0..cnvtxs {
                    bestwhere[k as usize] = cg.where_[k as usize];
                }
            }
        }

        // HiGHS: ometis.c lines 343-344: if (mincut == 0) break;
        if mincut == 0 {
            break;
        }

        // HiGHS: ometis.c lines 346-347: if (i < nruns-1) FreeRData(cgraph);
        if i < nruns - 1 {
            FreeRData(pool.get_mut(cgraph));
        }

        i += 1;
    }

    // HiGHS: ometis.c lines 350-351: if (mincut != cgraph->mincut) icopy(cgraph->nvtxs, bestwhere, cgraph->where);
    if mincut != pool.get(cgraph).mincut {
        let cg = pool.get_mut(cgraph);
        for k in 0..cnvtxs {
            cg.where_[k as usize] = bestwhere[k as usize];
        }
    }

    // HiGHS: ometis.c line 353: WCOREPOP;
    wspace::wspacepop(ctrl);

    // HiGHS: ometis.c line 355: Refine2WayNode(ctrl, graph, cgraph);
    Refine2WayNode(ctrl, pool, graph, cgraph);
}

/// HiGHS: ometis.c::MlevelNodeBisectionL1 (lines 363-380)
// PORT NOTE: the C if/else-if bound on CoarsenTo (lines 368-371) is kept as written.
#[allow(clippy::manual_clamp)]
pub fn MlevelNodeBisectionL1(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    niparts: idx_t,
) {
    // HiGHS: ometis.c line 367: ctrl->CoarsenTo = graph->nvtxs/8;
    ctrl.CoarsenTo = pool.get(graph).nvtxs / 8;
    // HiGHS: ometis.c lines 368-371
    if ctrl.CoarsenTo > 100 {
        ctrl.CoarsenTo = 100;
    } else if ctrl.CoarsenTo < 40 {
        ctrl.CoarsenTo = 40;
    }

    // HiGHS: ometis.c line 373: cgraph = CoarsenGraph(ctrl, graph);
    let cgraph = CoarsenGraph(ctrl, pool, graph);

    // HiGHS: ometis.c line 375: niparts = gk_max(1, (cgraph->nvtxs <= ctrl->CoarsenTo ? niparts/2: niparts));
    // gk_max is a C ternary: if a > b { a } else { b }
    {
        let a: idx_t = 1;
        let b: idx_t = if pool.get(cgraph).nvtxs <= ctrl.CoarsenTo {
            niparts / 2
        } else {
            niparts
        };
        let niparts = if a > b { a } else { b };

        // HiGHS: ometis.c line 377: InitSeparator(ctrl, cgraph, niparts);
        {
            let cg = pool.get_mut(cgraph);
            InitSeparator(ctrl, cg, niparts);
        }

        // HiGHS: ometis.c line 379: Refine2WayNode(ctrl, graph, cgraph);
        Refine2WayNode(ctrl, pool, graph, cgraph);
    }
}

/// HiGHS: ometis.c::SplitGraphOrder (lines 390-494)
pub fn SplitGraphOrder(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    r_lgraph: &mut Option<GraphId>,
    r_rgraph: &mut Option<GraphId>,
) {
    // HiGHS: ometis.c line 401: WCOREPUSH;
    wspace::wspacepush(ctrl);

    // HiGHS: ometis.c lines 403-411: local aliases
    let nvtxs = pool.get(graph).nvtxs;

    // HiGHS: ometis.c line 413: rename = iwspacemalloc(ctrl, nvtxs);
    let mut rename = wspace::iwspacemalloc(ctrl, nvtxs);

    // HiGHS: ometis.c line 415: snvtxs[0] = snvtxs[1] = snvtxs[2] = snedges[0] = snedges[1] = snedges[2] = 0;
    let mut snvtxs = [0 as idx_t; 3];
    let mut snedges = [0 as idx_t; 3];

    // HiGHS: ometis.c lines 416-420
    {
        let g = pool.get(graph);
        for i in 0..nvtxs {
            let k = g.where_[i as usize];
            rename[i as usize] = snvtxs[k as usize];
            snvtxs[k as usize] += 1;
            snedges[k as usize] += g.xadj_eff()[(i + 1) as usize] - g.xadj_eff()[i as usize];
        }
    }

    // HiGHS: ometis.c lines 422-427: lgraph = SetupSplitGraph(graph, snvtxs[0], snedges[0]);
    let lgraph = {
        let g = pool.get(graph);
        SetupSplitGraph(g, snvtxs[0], snedges[0])
    };
    let lgraph_id = pool.alloc(lgraph);

    // HiGHS: ometis.c lines 429-434: rgraph = SetupSplitGraph(graph, snvtxs[1], snedges[1]);
    let rgraph = {
        let g = pool.get(graph);
        SetupSplitGraph(g, snvtxs[1], snedges[1])
    };
    let rgraph_id = pool.alloc(rgraph);

    // HiGHS: ometis.c lines 436-441: Go and use bndptr to also mark the boundary nodes in the two partitions
    {
        let nbnd = pool.get(graph).nbnd;
        for ii in 0..nbnd {
            let i = pool.get(graph).bndind[ii as usize];
            let istart = pool.get(graph).xadj_eff()[i as usize];
            let iend = pool.get(graph).xadj_eff()[(i + 1) as usize];
            for j in istart..iend {
                let adj = pool.get(graph).adjncy_eff()[j as usize];
                pool.get_mut(graph).bndptr[adj as usize] = 1;
            }
        }
    }

    // HiGHS: ometis.c line 443: snvtxs[0] = snvtxs[1] = snedges[0] = snedges[1] = 0;
    snvtxs[0] = 0;
    snvtxs[1] = 0;
    snedges[0] = 0;
    snedges[1] = 0;

    // HiGHS: ometis.c line 444: sxadj[0][0] = sxadj[1][0] = 0;
    pool.get_mut(lgraph_id).xadj[0] = 0;
    pool.get_mut(rgraph_id).xadj[0] = 0;

    // HiGHS: ometis.c lines 445-471
    for i in 0..nvtxs {
        let mypart = pool.get(graph).where_[i as usize];
        if mypart == 2 {
            continue;
        }

        let istart = pool.get(graph).xadj_eff()[i as usize];
        let iend = pool.get(graph).xadj_eff()[(i + 1) as usize];

        if pool.get(graph).bndptr[i as usize] == -1 {
            // HiGHS: ometis.c lines 452-455: This is an interior vertex
            // auxadjncy = sadjncy[mypart] + snedges[mypart] - istart;
            // for(j=istart; j<iend; j++) auxadjncy[j] = adjncy[j];
            // snedges[mypart] += iend-istart;
            let offset = snedges[mypart as usize] - istart;
            {
                // Copy the needed adjncy values out of the immutable borrow first.
                let adjncy_vals: Vec<idx_t> = {
                    let g = pool.get(graph);
                    (istart..iend).map(|j| g.adjncy_eff()[j as usize]).collect()
                };
                let lg_or_rg = if mypart == 0 {
                    pool.get_mut(lgraph_id)
                } else {
                    pool.get_mut(rgraph_id)
                };
                for (idx, j) in (istart..iend).enumerate() {
                    lg_or_rg.adjncy[(offset + j) as usize] = adjncy_vals[idx];
                }
            }
            snedges[mypart as usize] += iend - istart;
        } else {
            // HiGHS: ometis.c lines 458-466: boundary vertex
            // auxadjncy = sadjncy[mypart];
            // l = snedges[mypart];
            // for (j=istart; j<iend; j++) { k = adjncy[j]; if (where[k] == mypart) auxadjncy[l++] = k; }
            // snedges[mypart] = l;
            let mut l = snedges[mypart as usize];
            {
                // Copy the needed adjncy and where values out of the immutable borrow first.
                let adjncy_vals: Vec<idx_t> = {
                    let g = pool.get(graph);
                    (istart..iend).map(|j| g.adjncy_eff()[j as usize]).collect()
                };
                let where_vals: Vec<idx_t> = {
                    let g = pool.get(graph);
                    (istart..iend)
                        .map(|j| g.where_[g.adjncy_eff()[j as usize] as usize])
                        .collect()
                };
                let lg_or_rg = if mypart == 0 {
                    pool.get_mut(lgraph_id)
                } else {
                    pool.get_mut(rgraph_id)
                };
                for (idx, _j) in (istart..iend).enumerate() {
                    let k = adjncy_vals[idx];
                    if where_vals[idx] == mypart {
                        lg_or_rg.adjncy[l as usize] = k;
                        l += 1;
                    }
                }
            }
            snedges[mypart as usize] = l;
        }

        // HiGHS: ometis.c lines 468-470
        {
            // Copy the needed vwgt and label values out of the immutable borrow first.
            let vwgt_val = pool.get(graph).vwgt[i as usize];
            let label_val = pool.get(graph).label[i as usize];
            let snvtxs_mypart = snvtxs[mypart as usize];
            let lg_or_rg = if mypart == 0 {
                pool.get_mut(lgraph_id)
            } else {
                pool.get_mut(rgraph_id)
            };
            lg_or_rg.vwgt[snvtxs_mypart as usize] = vwgt_val;
            lg_or_rg.label[snvtxs_mypart as usize] = label_val;
            snvtxs[mypart as usize] += 1;
            lg_or_rg.xadj[snvtxs[mypart as usize] as usize] = snedges[mypart as usize];
        }
    }

    // HiGHS: ometis.c lines 473-480
    for mypart in 0..2 {
        let iend = snedges[mypart];
        {
            let lg_or_rg = if mypart == 0 {
                pool.get_mut(lgraph_id)
            } else {
                pool.get_mut(rgraph_id)
            };
            // iset(iend, 1, sadjwgt[mypart]);
            for k in 0..iend {
                lg_or_rg.adjwgt[k as usize] = 1;
            }
        }
        {
            let lg_or_rg = if mypart == 0 {
                pool.get_mut(lgraph_id)
            } else {
                pool.get_mut(rgraph_id)
            };
            // auxadjncy = sadjncy[mypart];
            // for (i=0; i<iend; i++) auxadjncy[i] = rename[auxadjncy[i]];
            for i in 0..iend {
                lg_or_rg.adjncy[i as usize] = rename[lg_or_rg.adjncy[i as usize] as usize];
            }
        }
    }

    // HiGHS: ometis.c lines 482-485
    pool.get_mut(lgraph_id).nvtxs = snvtxs[0];
    pool.get_mut(lgraph_id).nedges = snedges[0];
    pool.get_mut(rgraph_id).nvtxs = snvtxs[1];
    pool.get_mut(rgraph_id).nedges = snedges[1];

    // HiGHS: ometis.c lines 487-488
    SetupGraph_tvwgt(pool.get_mut(lgraph_id));
    SetupGraph_tvwgt(pool.get_mut(rgraph_id));

    // HiGHS: ometis.c lines 490-491
    *r_lgraph = Some(lgraph_id);
    *r_rgraph = Some(rgraph_id);

    // HiGHS: ometis.c line 493: WCOREPOP;
    wspace::wspacepop(ctrl);
}

/// HiGHS: ometis.c::SplitGraphOrderCC (lines 516-608)
pub fn SplitGraphOrderCC(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    ncmps: idx_t,
    cptr: &[idx_t],
    cind: &[idx_t],
) -> Vec<GraphId> {
    // HiGHS: ometis.c line 527: WCOREPUSH;
    wspace::wspacepush(ctrl);

    // HiGHS: ometis.c lines 529-537: local aliases
    let nvtxs = pool.get(graph).nvtxs;

    // HiGHS: ometis.c lines 539-544: Go and use bndptr to also mark the boundary nodes in the two partitions
    {
        let nbnd = pool.get(graph).nbnd;
        for ii in 0..nbnd {
            let i = pool.get(graph).bndind[ii as usize];
            let istart = pool.get(graph).xadj_eff()[i as usize];
            let iend = pool.get(graph).xadj_eff()[(i + 1) as usize];
            for j in istart..iend {
                let adj = pool.get(graph).adjncy_eff()[j as usize];
                pool.get_mut(graph).bndptr[adj as usize] = 1;
            }
        }
    }

    // HiGHS: ometis.c line 546: rename = iwspacemalloc(ctrl, nvtxs);
    let mut rename = wspace::iwspacemalloc(ctrl, nvtxs);

    // HiGHS: ometis.c line 548: sgraphs = (graph_t **)malloc(sizeof(graph_t *)*ncmps);
    let mut sgraphs: Vec<GraphId> = Vec::with_capacity(ncmps as usize);

    // The C code mutates cind in place via irandArrayPermute. With &[idx_t] we cannot,
    // so we make a full mutable copy and use it throughout.
    let mut cind_mut = cind.to_vec();

    // HiGHS: ometis.c lines 551-603: Go and split the graph a component at a time
    for iii in 0..ncmps {
        // HiGHS: ometis.c line 552: irandArrayPermute(cptr[iii+1]-cptr[iii], cind+cptr[iii], cptr[iii+1]-cptr[iii], 0, &ctrl->rng_state);
        {
            let n = cptr[(iii + 1) as usize] - cptr[iii as usize];
            let offset = cptr[iii as usize] as usize;
            let len = n as usize;
            super::super::gklib::irandArrayPermute(
                n,
                &mut cind_mut[offset..offset + len],
                n,
                0,
                &mut ctrl.rng_state,
            );
        }

        // HiGHS: ometis.c line 553: snvtxs = snedges = 0;
        let mut snvtxs: idx_t = 0;
        let mut snedges: idx_t = 0;

        // HiGHS: ometis.c lines 554-558
        {
            let g = pool.get(graph);
            for j in cptr[iii as usize]..cptr[(iii + 1) as usize] {
                let i = cind_mut[j as usize];
                rename[i as usize] = snvtxs;
                snvtxs += 1;
                snedges += g.xadj_eff()[(i + 1) as usize] - g.xadj_eff()[i as usize];
            }
        }

        // HiGHS: ometis.c line 560: sgraphs[iii] = SetupSplitGraph(graph, snvtxs, snedges);
        let sgraph = {
            let g = pool.get(graph);
            SetupSplitGraph(g, snvtxs, snedges)
        };
        let sgraph_id = pool.alloc(sgraph);
        sgraphs.push(sgraph_id);

        // HiGHS: ometis.c line 568: snvtxs = snedges = sxadj[0] = 0;
        snvtxs = 0;
        snedges = 0;
        pool.get_mut(sgraph_id).xadj[0] = 0;

        // HiGHS: ometis.c lines 569-593
        for ii in cptr[iii as usize]..cptr[(iii + 1) as usize] {
            let i = cind_mut[ii as usize];

            let istart = pool.get(graph).xadj_eff()[i as usize];
            let iend = pool.get(graph).xadj_eff()[(i + 1) as usize];

            if pool.get(graph).bndptr[i as usize] == -1 {
                // HiGHS: ometis.c lines 574-578: interior vertex
                // auxadjncy = sadjncy + snedges - istart;
                // for(j=istart; j<iend; j++) auxadjncy[j] = adjncy[j];
                // snedges += iend-istart;
                let offset = snedges - istart;
                {
                    // Copy the needed adjncy values out of the immutable borrow first.
                    let adjncy_vals: Vec<idx_t> = {
                        let g = pool.get(graph);
                        (istart..iend).map(|j| g.adjncy_eff()[j as usize]).collect()
                    };
                    let sg = pool.get_mut(sgraph_id);
                    for (idx, j) in (istart..iend).enumerate() {
                        sg.adjncy[(offset + j) as usize] = adjncy_vals[idx];
                    }
                }
                snedges += iend - istart;
            } else {
                // HiGHS: ometis.c lines 580-588: boundary vertex
                // l = snedges;
                // for (j=istart; j<iend; j++) { k = adjncy[j]; if (where[k] != 2) sadjncy[l++] = k; }
                // snedges = l;
                let mut l = snedges;
                {
                    // Copy the needed adjncy and where values out of the immutable borrow first.
                    let adjncy_vals: Vec<idx_t> = {
                        let g = pool.get(graph);
                        (istart..iend).map(|j| g.adjncy_eff()[j as usize]).collect()
                    };
                    let where_vals: Vec<idx_t> = {
                        let g = pool.get(graph);
                        (istart..iend)
                            .map(|j| g.where_[g.adjncy_eff()[j as usize] as usize])
                            .collect()
                    };
                    let sg = pool.get_mut(sgraph_id);
                    for (idx, _j) in (istart..iend).enumerate() {
                        let k = adjncy_vals[idx];
                        if where_vals[idx] != 2 {
                            sg.adjncy[l as usize] = k;
                            l += 1;
                        }
                    }
                }
                snedges = l;
            }

            // HiGHS: ometis.c lines 590-592
            {
                // Copy the needed vwgt and label values out of the immutable borrow first.
                let vwgt_val = pool.get(graph).vwgt[i as usize];
                let label_val = pool.get(graph).label[i as usize];
                let sg = pool.get_mut(sgraph_id);
                sg.vwgt[snvtxs as usize] = vwgt_val;
                sg.label[snvtxs as usize] = label_val;
                snvtxs += 1;
                sg.xadj[snvtxs as usize] = snedges;
            }
        }

        // HiGHS: ometis.c line 595: iset(snedges, 1, sadjwgt);
        {
            let sg = pool.get_mut(sgraph_id);
            for k in 0..snedges {
                sg.adjwgt[k as usize] = 1;
            }
        }

        // HiGHS: ometis.c lines 596-597: for (i=0; i<snedges; i++) sadjncy[i] = rename[sadjncy[i]];
        {
            let sg = pool.get_mut(sgraph_id);
            for i in 0..snedges {
                sg.adjncy[i as usize] = rename[sg.adjncy[i as usize] as usize];
            }
        }

        // HiGHS: ometis.c lines 599-600
        pool.get_mut(sgraph_id).nvtxs = snvtxs;
        pool.get_mut(sgraph_id).nedges = snedges;

        // HiGHS: ometis.c line 602: SetupGraph_tvwgt(sgraphs[iii]);
        SetupGraph_tvwgt(pool.get_mut(sgraph_id));
    }

    // HiGHS: ometis.c line 605: WCOREPOP;
    wspace::wspacepop(ctrl);

    // HiGHS: ometis.c line 607: return sgraphs;
    sgraphs
}

/// HiGHS: ometis.c::MMDOrder (lines 615-656)
/// Shifts `graph->xadj` / `adjncy` to 1-based in place, runs `genmmd`,
/// shifts back.
pub fn MMDOrder(ctrl: &mut ctrl_t, graph: &mut graph_t, order: &mut [idx_t], lastvtx: idx_t) {
    // HiGHS: ometis.c line 621: WCOREPUSH;
    wspace::wspacepush(ctrl);

    // HiGHS: ometis.c lines 623-625
    let nvtxs = graph.nvtxs;

    // HiGHS: ometis.c lines 627-632: Relabel the vertices so that it starts from 1
    let k = graph.xadj[nvtxs as usize];
    for i in 0..k {
        graph.adjncy[i as usize] += 1;
    }
    for i in 0..(nvtxs + 1) {
        graph.xadj[i as usize] += 1;
    }

    // HiGHS: ometis.c lines 634-639: workspace allocations
    let mut perm = wspace::iwspacemalloc(ctrl, nvtxs + 5);
    let mut iperm = wspace::iwspacemalloc(ctrl, nvtxs + 5);
    let mut head = wspace::iwspacemalloc(ctrl, nvtxs + 5);
    let mut qsize = wspace::iwspacemalloc(ctrl, nvtxs + 5);
    let mut list = wspace::iwspacemalloc(ctrl, nvtxs + 5);
    let mut marker = wspace::iwspacemalloc(ctrl, nvtxs + 5);

    // HiGHS: ometis.c line 641: genmmd(nvtxs, xadj, adjncy, iperm, perm, 1, head, qsize, list, marker, IDX_MAX, &nofsub);
    let mut nofsub: idx_t = 0;
    genmmd(
        nvtxs,
        &mut graph.xadj,
        &mut graph.adjncy,
        &mut iperm,
        &mut perm,
        1,
        &mut head,
        &mut qsize,
        &mut list,
        &mut marker,
        idx_t::MAX,
        &mut nofsub,
    );

    // HiGHS: ometis.c lines 643-646
    let firstvtx = lastvtx - nvtxs;
    for i in 0..nvtxs {
        order[graph.label[i as usize] as usize] = firstvtx + iperm[i as usize] - 1;
    }

    // HiGHS: ometis.c lines 648-653: Relabel the vertices so that it starts from 0
    for i in 0..(nvtxs + 1) {
        graph.xadj[i as usize] -= 1;
    }
    let k = graph.xadj[nvtxs as usize];
    for i in 0..k {
        graph.adjncy[i as usize] -= 1;
    }

    // HiGHS: ometis.c line 655: WCOREPOP;
    wspace::wspacepop(ctrl);
}
