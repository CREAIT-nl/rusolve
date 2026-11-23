//! HiGHS source: extern/metis/libmetis/coarsen.c (lines 1-415)
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Multilevel coarsening drivers and the RM / SHEM matchings. The random
//! visit order comes from `irandArrayPermute` / `BucketSortKeysInc` on
//! `ctrl->rng_state` and must match the C call for call.

pub mod part2;

pub use part2::*;

use super::{GraphId, GraphPool, ctrl_t, graph_t, idx_t};

use super::gklib::{gk_errexit, imalloc, irandArrayPermute, iset};
use super::mcutil::{BetterVBalance, ivecaxpylez, ivecle};
use super::refine::BucketSortKeysInc;
use super::wspace::{iwspacemalloc, wspacepop, wspacepush};
use super::{COARSEN_FRACTION, METIS_CTYPE_RM, METIS_CTYPE_SHEM, UNMATCHED, UNMATCHEDFOR2HOP};

/// HiGHS: coarsen.c::CoarsenGraph (lines 22-74)
/// Returns the coarsest graph of the chain built below `graph`.
pub fn CoarsenGraph(ctrl: &mut ctrl_t, pool: &mut GraphPool, graph: GraphId) -> GraphId {
    // HiGHS: coarsen.c lines 22-74
    let mut graph = graph;
    let mut eqewgts: idx_t;
    // HiGHS: coarsen.c line 24 `level` is written at line 64 and never read; not kept

    // determine if the weights on the edges are all the same
    // HiGHS: coarsen.c lines 28-33
    eqewgts = 1;
    {
        let g = pool.get(graph);
        let adjwgt = &g.adjwgt;
        let nedges = g.nedges;
        let mut i: idx_t = 1;
        while i < nedges {
            if adjwgt[0] != adjwgt[i as usize] {
                eqewgts = 0;
                break;
            }
            i += 1;
        }
    }

    // set the maximum allowed coarsest vertex weight
    // HiGHS: coarsen.c lines 36-37
    {
        let g = pool.get(graph);
        let ncon = g.ncon;
        let tvwgt = &g.tvwgt;
        for i in 0..ncon {
            ctrl.maxvwgt[i as usize] =
                ((1.5 * tvwgt[i as usize] as f64) / ctrl.CoarsenTo as f64) as idx_t;
        }
    }

    // HiGHS: coarsen.c lines 39-69 (do-while)
    loop {
        // IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));
        // lines 40 dropped (dbglvl off in oracle build)

        // allocate memory for cmap, if it has not already been done due to multiple cuts
        // HiGHS: coarsen.c lines 44-45
        if pool.get(graph).cmap.is_empty() {
            let nvtxs = pool.get(graph).nvtxs;
            pool.get_mut(graph).cmap = imalloc(nvtxs as usize);
        }

        // determine which matching scheme you will use
        // HiGHS: coarsen.c lines 48-60
        match ctrl.ctype {
            METIS_CTYPE_RM => {
                Match_RM(ctrl, pool, graph);
            }
            METIS_CTYPE_SHEM => {
                let g = pool.get(graph);
                let nedges = g.nedges;
                if eqewgts != 0 || nedges == 0 {
                    Match_RM(ctrl, pool, graph);
                } else {
                    Match_SHEM(ctrl, pool, graph);
                }
            }
            _ => {
                gk_errexit(&format!("Unknown ctype: {}\n", ctrl.ctype));
            }
        }

        // HiGHS: coarsen.c line 62
        graph = pool.get(graph).coarser.unwrap();
        // HiGHS: coarsen.c line 63
        eqewgts = 0;
        // HiGHS: coarsen.c line 64 (`level++`, never read)

        // HiGHS: coarsen.c lines 67-69 (while condition)
        let g = pool.get(graph);
        let nvtxs = g.nvtxs;
        let nedges = g.nedges;
        let finer_nvtxs = pool.get(g.finer.unwrap()).nvtxs;

        if !(nvtxs > ctrl.CoarsenTo
            && (nvtxs as f64) < COARSEN_FRACTION * (finer_nvtxs as f64)
            && nedges > nvtxs / 2)
        {
            break;
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));
    // line 71 dropped (dbglvl off in oracle build)

    graph
}

/// HiGHS: coarsen.c::CoarsenGraphNlevels (lines 82-135)
pub fn CoarsenGraphNlevels(
    ctrl: &mut ctrl_t,
    pool: &mut GraphPool,
    graph: GraphId,
    nlevels: idx_t,
) -> GraphId {
    // HiGHS: coarsen.c lines 82-135
    let mut graph = graph;
    let mut eqewgts: idx_t;

    // determine if the weights on the edges are all the same
    // HiGHS: coarsen.c lines 88-93
    eqewgts = 1;
    {
        let g = pool.get(graph);
        let adjwgt = &g.adjwgt;
        let nedges = g.nedges;
        let mut i: idx_t = 1;
        while i < nedges {
            if adjwgt[0] != adjwgt[i as usize] {
                eqewgts = 0;
                break;
            }
            i += 1;
        }
    }

    // set the maximum allowed coarsest vertex weight
    // HiGHS: coarsen.c lines 96-97
    {
        let g = pool.get(graph);
        let ncon = g.ncon;
        let tvwgt = &g.tvwgt;
        for i in 0..ncon {
            ctrl.maxvwgt[i as usize] =
                ((1.5 * tvwgt[i as usize] as f64) / ctrl.CoarsenTo as f64) as idx_t;
        }
    }

    // HiGHS: coarsen.c lines 99-130 (for loop)
    let mut level: idx_t = 0;
    while level < nlevels {
        // IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));
        // line 100 dropped (dbglvl off in oracle build)

        // allocate memory for cmap, if it has not already been done due to multiple cuts
        // HiGHS: coarsen.c lines 104-105
        if pool.get(graph).cmap.is_empty() {
            let nvtxs = pool.get(graph).nvtxs;
            pool.get_mut(graph).cmap = imalloc(nvtxs as usize);
        }

        // determine which matching scheme you will use
        // HiGHS: coarsen.c lines 108-120
        match ctrl.ctype {
            METIS_CTYPE_RM => {
                Match_RM(ctrl, pool, graph);
            }
            METIS_CTYPE_SHEM => {
                let g = pool.get(graph);
                let nedges = g.nedges;
                if eqewgts != 0 || nedges == 0 {
                    Match_RM(ctrl, pool, graph);
                } else {
                    Match_SHEM(ctrl, pool, graph);
                }
            }
            _ => {
                gk_errexit(&format!("Unknown ctype: {}\n", ctrl.ctype));
            }
        }

        // HiGHS: coarsen.c line 122
        graph = pool.get(graph).coarser.unwrap();
        // HiGHS: coarsen.c line 123
        eqewgts = 0;

        // HiGHS: coarsen.c lines 126-129 (break condition)
        let g = pool.get(graph);
        let nvtxs = g.nvtxs;
        let nedges = g.nedges;
        let finer_nvtxs = pool.get(g.finer.unwrap()).nvtxs;

        if nvtxs < ctrl.CoarsenTo
            || (nvtxs as f64) > COARSEN_FRACTION * (finer_nvtxs as f64)
            || nedges < nvtxs / 2
        {
            break;
        }

        level += 1;
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_COARSEN, PrintCGraphStats(ctrl, graph));
    // line 132 dropped (dbglvl off in oracle build)

    graph
}

/// HiGHS: coarsen.c::Match_RM (lines 143-272)
/// Ends with `CreateCoarseGraph`, which links a new coarser graph.
pub fn Match_RM(ctrl: &mut ctrl_t, pool: &mut GraphPool, graph: GraphId) -> idx_t {
    // HiGHS: coarsen.c lines 143-272
    let mut i: idx_t;
    let mut pi: idx_t;
    let mut j: idx_t;
    let mut k: idx_t;
    let mut nvtxs: idx_t;
    let mut ncon: idx_t;
    let mut cnvtxs: idx_t;
    let mut maxidx: idx_t;
    let mut last_unmatched: idx_t;
    let mut avgdegree: idx_t;
    let mut bnum: idx_t;
    let mut nunmatched: usize = 0;

    // WCOREPUSH
    // HiGHS: coarsen.c line 152
    wspacepush(ctrl);

    // HiGHS: coarsen.c lines 155-161
    nvtxs = pool.get(graph).nvtxs;
    ncon = pool.get(graph).ncon;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
    // vwgt = graph->vwgt;
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
    // adjwgt = graph->adjwgt;
    // cmap = graph->cmap;

    // HiGHS: coarsen.c line 163
    // maxvwgt = ctrl->maxvwgt;

    // HiGHS: coarsen.c lines 165-168
    let mut match_: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    iset(nvtxs as usize, UNMATCHED, &mut match_);
    let mut perm: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    let mut tperm: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    let mut degrees: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);

    // Determine a "random" traversal order that is biased towards low-degree vertices
    // HiGHS: coarsen.c line 172
    irandArrayPermute(nvtxs, &mut tperm, nvtxs / 8, 1, &mut ctrl.rng_state);

    // HiGHS: coarsen.c line 174
    {
        let xadj = pool.get(graph).xadj_eff();
        avgdegree = (4.0 * ((xadj[nvtxs as usize] / nvtxs) as f64)) as idx_t;
    }

    // HiGHS: coarsen.c lines 175-178
    {
        let xadj = pool.get(graph).xadj_eff();
        for i in 0..nvtxs {
            bnum = ((1 + xadj[(i + 1) as usize] - xadj[i as usize]) as f64).sqrt() as idx_t;
            degrees[i as usize] = if bnum > avgdegree { avgdegree } else { bnum };
        }
    }

    // HiGHS: coarsen.c line 179
    BucketSortKeysInc(ctrl, nvtxs, avgdegree, &degrees, &tperm, &mut perm);

    // Traverse the vertices and compute the matching
    // HiGHS: coarsen.c lines 183-246
    cnvtxs = 0;
    last_unmatched = 0;
    {
        // Take cmap out for the traversal
        let mut cmap = std::mem::take(&mut pool.get_mut(graph).cmap);

        for pi in 0..nvtxs {
            i = perm[pi as usize];

            if match_[i as usize] == UNMATCHED {
                maxidx = i;

                // HiGHS: coarsen.c line 189
                let cond = if ncon == 1 {
                    pool.get(graph).vwgt[i as usize] < ctrl.maxvwgt[0]
                } else {
                    let vwgt = &pool.get(graph).vwgt;
                    ivecle(ncon, &vwgt[(i * ncon) as usize..], &ctrl.maxvwgt) != 0
                };

                if cond {
                    // Deal with island vertices. Find a non-island and match it with.
                    // The matching ignores ctrl->maxvwgt requirements
                    // HiGHS: coarsen.c line 192
                    let xadj = pool.get(graph).xadj_eff();
                    if xadj[i as usize] == xadj[(i + 1) as usize] {
                        // HiGHS: coarsen.c line 193
                        last_unmatched = (if pi > last_unmatched {
                            pi
                        } else {
                            last_unmatched
                        }) + 1;
                        // HiGHS: coarsen.c lines 194-200
                        while last_unmatched < nvtxs {
                            j = perm[last_unmatched as usize];
                            if match_[j as usize] == UNMATCHED {
                                maxidx = j;
                                break;
                            }
                            last_unmatched += 1;
                        }
                    } else {
                        // Find a random matching, subject to maxvwgt constraints
                        // HiGHS: coarsen.c lines 204-236
                        if ncon == 1 {
                            // single constraint version
                            // HiGHS: coarsen.c lines 206-212
                            let xadj = pool.get(graph).xadj_eff();
                            let adjncy = pool.get(graph).adjncy_eff();
                            let vwgt = &pool.get(graph).vwgt;
                            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                                k = adjncy[j as usize];
                                if match_[k as usize] == UNMATCHED
                                    && vwgt[i as usize] + vwgt[k as usize] <= ctrl.maxvwgt[0]
                                {
                                    maxidx = k;
                                    break;
                                }
                            }

                            // If it did not match, record for a 2-hop matching.
                            // HiGHS: coarsen.c lines 215-218
                            if maxidx == i && 2 * vwgt[i as usize] < ctrl.maxvwgt[0] {
                                nunmatched += 1;
                                maxidx = UNMATCHED;
                            }
                        } else {
                            // multi-constraint version
                            // HiGHS: coarsen.c lines 222-229
                            let xadj = pool.get(graph).xadj_eff();
                            let adjncy = pool.get(graph).adjncy_eff();
                            let vwgt = &pool.get(graph).vwgt;
                            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                                k = adjncy[j as usize];
                                if match_[k as usize] == UNMATCHED
                                    && ivecaxpylez(
                                        ncon,
                                        1,
                                        &vwgt[(i * ncon) as usize..],
                                        &vwgt[(k * ncon) as usize..],
                                        &ctrl.maxvwgt,
                                    ) != 0
                                {
                                    maxidx = k;
                                    break;
                                }
                            }

                            // If it did not match, record for a 2-hop matching.
                            // HiGHS: coarsen.c lines 232-235
                            if maxidx == i
                                && ivecaxpylez(
                                    ncon,
                                    2,
                                    &vwgt[(i * ncon) as usize..],
                                    &vwgt[(i * ncon) as usize..],
                                    &ctrl.maxvwgt,
                                ) != 0
                            {
                                nunmatched += 1;
                                maxidx = UNMATCHED;
                            }
                        }
                    }
                }

                // HiGHS: coarsen.c lines 240-244
                if maxidx != UNMATCHED {
                    let c = cnvtxs;
                    cmap[i as usize] = c;
                    cmap[maxidx as usize] = c;
                    cnvtxs += 1;
                    match_[i as usize] = maxidx;
                    match_[maxidx as usize] = i;
                }
            }
        }

        // Put cmap back
        pool.get_mut(graph).cmap = cmap;
    }

    // see if a 2-hop matching is required/allowed
    // HiGHS: coarsen.c lines 249-250
    if ctrl.no2hop == 0 && (nunmatched as f64) > UNMATCHEDFOR2HOP * (nvtxs as f64) {
        // HiGHS: coarsen.c line 250 / 393; the result is overwritten by the loop below
        Match_2Hop(
            ctrl,
            pool.get_mut(graph),
            &mut perm,
            &mut match_,
            cnvtxs,
            nunmatched,
        );
    }

    // match the final unmatched vertices with themselves and reorder the vertices
    // of the coarse graph for memory-friendly contraction
    // HiGHS: coarsen.c lines 255-264
    {
        let mut cmap = std::mem::take(&mut pool.get_mut(graph).cmap);
        cnvtxs = 0;
        for i in 0..nvtxs {
            if match_[i as usize] == UNMATCHED {
                match_[i as usize] = i;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
            } else {
                if i <= match_[i as usize] {
                    let c = cnvtxs;
                    cmap[i as usize] = c;
                    cmap[match_[i as usize] as usize] = c;
                    cnvtxs += 1;
                }
            }
        }
        pool.get_mut(graph).cmap = cmap;
    }

    // HiGHS: coarsen.c line 267
    CreateCoarseGraph(ctrl, pool, graph, cnvtxs, &match_);

    // WCOREPOP
    // HiGHS: coarsen.c line 269
    wspacepop(ctrl);

    cnvtxs
}

/// HiGHS: coarsen.c::Match_SHEM (lines 281-415)
pub fn Match_SHEM(ctrl: &mut ctrl_t, pool: &mut GraphPool, graph: GraphId) -> idx_t {
    // HiGHS: coarsen.c lines 281-415
    let mut i: idx_t;
    let mut pi: idx_t;
    let mut j: idx_t;
    let mut k: idx_t;
    let mut nvtxs: idx_t;
    let mut ncon: idx_t;
    let mut cnvtxs: idx_t;
    let mut maxidx: idx_t;
    let mut maxwgt: idx_t;
    let mut last_unmatched: idx_t;
    let mut avgdegree: idx_t;
    let mut bnum: idx_t;
    let mut nunmatched: usize = 0;

    // WCOREPUSH
    // HiGHS: coarsen.c line 290
    wspacepush(ctrl);

    // HiGHS: coarsen.c lines 293-299
    nvtxs = pool.get(graph).nvtxs;
    ncon = pool.get(graph).ncon;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj;
    // vwgt = graph->vwgt;
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy;
    // adjwgt = graph->adjwgt;
    // cmap = graph->cmap;

    // HiGHS: coarsen.c line 301
    // maxvwgt = ctrl->maxvwgt;

    // HiGHS: coarsen.c lines 303-306
    let mut match_: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    iset(nvtxs as usize, UNMATCHED, &mut match_);
    let mut perm: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    let mut tperm: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);
    let mut degrees: Vec<idx_t> = iwspacemalloc(ctrl, nvtxs);

    // Determine a "random" traversal order that is biased towards low-degree vertices
    // HiGHS: coarsen.c line 309
    irandArrayPermute(nvtxs, &mut tperm, nvtxs / 8, 1, &mut ctrl.rng_state);

    // HiGHS: coarsen.c line 311
    {
        let xadj = pool.get(graph).xadj_eff();
        avgdegree = (4.0 * ((xadj[nvtxs as usize] / nvtxs) as f64)) as idx_t;
    }

    // HiGHS: coarsen.c lines 312-315
    {
        let xadj = pool.get(graph).xadj_eff();
        for i in 0..nvtxs {
            bnum = ((1 + xadj[(i + 1) as usize] - xadj[i as usize]) as f64).sqrt() as idx_t;
            degrees[i as usize] = if bnum > avgdegree { avgdegree } else { bnum };
        }
    }

    // HiGHS: coarsen.c line 316
    BucketSortKeysInc(ctrl, nvtxs, avgdegree, &degrees, &tperm, &mut perm);

    // Traverse the vertices and compute the matching
    // HiGHS: coarsen.c lines 320-389
    cnvtxs = 0;
    last_unmatched = 0;
    {
        // Take cmap out for the traversal
        let mut cmap = std::mem::take(&mut pool.get_mut(graph).cmap);

        for pi in 0..nvtxs {
            i = perm[pi as usize];

            if match_[i as usize] == UNMATCHED {
                maxidx = i;
                maxwgt = -1;

                // HiGHS: coarsen.c line 327
                let cond = if ncon == 1 {
                    pool.get(graph).vwgt[i as usize] < ctrl.maxvwgt[0]
                } else {
                    let vwgt = &pool.get(graph).vwgt;
                    ivecle(ncon, &vwgt[(i * ncon) as usize..], &ctrl.maxvwgt) != 0
                };

                if cond {
                    // Deal with island vertices. Find a non-island and match it with.
                    // The matching ignores ctrl->maxvwgt requirements
                    // HiGHS: coarsen.c line 330
                    let xadj = pool.get(graph).xadj_eff();
                    if xadj[i as usize] == xadj[(i + 1) as usize] {
                        // HiGHS: coarsen.c line 331
                        last_unmatched = (if pi > last_unmatched {
                            pi
                        } else {
                            last_unmatched
                        }) + 1;
                        // HiGHS: coarsen.c lines 332-338
                        while last_unmatched < nvtxs {
                            j = perm[last_unmatched as usize];
                            if match_[j as usize] == UNMATCHED {
                                maxidx = j;
                                break;
                            }
                            last_unmatched += 1;
                        }
                    } else {
                        // Find a heavy-edge matching, subject to maxvwgt constraints
                        // HiGHS: coarsen.c lines 342-379
                        if ncon == 1 {
                            // single constraint version
                            // HiGHS: coarsen.c lines 344-351
                            let xadj = pool.get(graph).xadj_eff();
                            let adjncy = pool.get(graph).adjncy_eff();
                            let adjwgt = &pool.get(graph).adjwgt;
                            let vwgt = &pool.get(graph).vwgt;
                            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                                k = adjncy[j as usize];
                                if match_[k as usize] == UNMATCHED
                                    && maxwgt < adjwgt[j as usize]
                                    && vwgt[i as usize] + vwgt[k as usize] <= ctrl.maxvwgt[0]
                                {
                                    maxidx = k;
                                    maxwgt = adjwgt[j as usize];
                                }
                            }

                            // If it did not match, record for a 2-hop matching.
                            // HiGHS: coarsen.c lines 354-357
                            if maxidx == i && 2 * vwgt[i as usize] < ctrl.maxvwgt[0] {
                                nunmatched += 1;
                                maxidx = UNMATCHED;
                            }
                        } else {
                            // multi-constraint version
                            // HiGHS: coarsen.c lines 361-372
                            let xadj = pool.get(graph).xadj_eff();
                            let adjncy = pool.get(graph).adjncy_eff();
                            let adjwgt = &pool.get(graph).adjwgt;
                            let vwgt = &pool.get(graph).vwgt;
                            let invtvwgt = &pool.get(graph).invtvwgt;
                            for j in xadj[i as usize]..xadj[(i + 1) as usize] {
                                k = adjncy[j as usize];
                                if match_[k as usize] == UNMATCHED
                                    && ivecaxpylez(
                                        ncon,
                                        1,
                                        &vwgt[(i * ncon) as usize..],
                                        &vwgt[(k * ncon) as usize..],
                                        &ctrl.maxvwgt,
                                    ) != 0
                                    && (maxwgt < adjwgt[j as usize]
                                        || (maxwgt == adjwgt[j as usize]
                                            && BetterVBalance(
                                                ncon,
                                                invtvwgt,
                                                &vwgt[(i * ncon) as usize..],
                                                &vwgt[(maxidx * ncon) as usize..],
                                                &vwgt[(k * ncon) as usize..],
                                            ) != 0))
                                {
                                    maxidx = k;
                                    maxwgt = adjwgt[j as usize];
                                }
                            }

                            // If it did not match, record for a 2-hop matching.
                            // HiGHS: coarsen.c lines 375-378
                            if maxidx == i
                                && ivecaxpylez(
                                    ncon,
                                    2,
                                    &vwgt[(i * ncon) as usize..],
                                    &vwgt[(i * ncon) as usize..],
                                    &ctrl.maxvwgt,
                                ) != 0
                            {
                                nunmatched += 1;
                                maxidx = UNMATCHED;
                            }
                        }
                    }
                }

                // HiGHS: coarsen.c lines 383-387
                if maxidx != UNMATCHED {
                    let c = cnvtxs;
                    cmap[i as usize] = c;
                    cmap[maxidx as usize] = c;
                    cnvtxs += 1;
                    match_[i as usize] = maxidx;
                    match_[maxidx as usize] = i;
                }
            }
        }

        // Put cmap back
        pool.get_mut(graph).cmap = cmap;
    }

    // see if a 2-hop matching is required/allowed
    // HiGHS: coarsen.c lines 392-393
    if ctrl.no2hop == 0 && (nunmatched as f64) > UNMATCHEDFOR2HOP * (nvtxs as f64) {
        // HiGHS: coarsen.c line 250 / 393; the result is overwritten by the loop below
        Match_2Hop(
            ctrl,
            pool.get_mut(graph),
            &mut perm,
            &mut match_,
            cnvtxs,
            nunmatched,
        );
    }

    // match the final unmatched vertices with themselves and reorder the vertices
    // of the coarse graph for memory-friendly contraction
    // HiGHS: coarsen.c lines 398-407
    {
        let mut cmap = std::mem::take(&mut pool.get_mut(graph).cmap);
        cnvtxs = 0;
        for i in 0..nvtxs {
            if match_[i as usize] == UNMATCHED {
                match_[i as usize] = i;
                cmap[i as usize] = cnvtxs;
                cnvtxs += 1;
            } else {
                if i <= match_[i as usize] {
                    let c = cnvtxs;
                    cmap[i as usize] = c;
                    cmap[match_[i as usize] as usize] = c;
                    cnvtxs += 1;
                }
            }
        }
        pool.get_mut(graph).cmap = cmap;
    }

    // HiGHS: coarsen.c line 410
    CreateCoarseGraph(ctrl, pool, graph, cnvtxs, &match_);

    // WCOREPOP
    // HiGHS: coarsen.c line 412
    wspacepop(ctrl);

    cnvtxs
}

#[cfg(test)]
mod tests {
    //! Unit tests. Hand-derived from coarsen.c
    //! lines 22-415, not oracle fixtures. Each graph is chosen so that the
    //! expected matching does not depend on the `irandArrayPermute` visit
    //! order, except `match_shem_ncon2_tie_uses_better_vbalance`, which picks a
    //! seed whose visit order starts at the vertex under test.
    use super::super::graph::{CreateGraph, SetupGraph_tvwgt};
    use super::*;

    fn graph_from(
        xadj: Vec<idx_t>,
        adjncy: Vec<idx_t>,
        adjwgt: Vec<idx_t>,
        ncon: idx_t,
        vwgt: Vec<idx_t>,
    ) -> graph_t {
        let mut g = CreateGraph();
        g.nvtxs = (xadj.len() - 1) as idx_t;
        g.nedges = adjncy.len() as idx_t;
        g.ncon = ncon;
        assert_eq!(vwgt.len(), (g.nvtxs * ncon) as usize);
        assert_eq!(adjwgt.len(), adjncy.len());
        g.xadj = xadj;
        g.adjncy = adjncy;
        g.adjwgt = adjwgt;
        g.vwgt = vwgt;
        SetupGraph_tvwgt(&mut g);
        g
    }

    fn ctrl_with(ncon: idx_t, maxvwgt: Vec<idx_t>) -> ctrl_t {
        let mut ctrl = ctrl_t::default();
        ctrl.ctype = METIS_CTYPE_RM;
        ctrl.CoarsenTo = 1;
        ctrl.maxvwgt = maxvwgt;
        assert_eq!(ctrl.maxvwgt.len(), ncon as usize);
        ctrl
    }

    /// Calls `matcher` on a pool holding `g` with `cmap` allocated (as
    /// `CoarsenGraph` lines 44-45 do) and returns (pool, fine id, cnvtxs).
    fn run_match(
        ctrl: &mut ctrl_t,
        mut g: graph_t,
        matcher: fn(&mut ctrl_t, &mut GraphPool, GraphId) -> idx_t,
    ) -> (GraphPool, GraphId, idx_t) {
        g.cmap = imalloc(g.nvtxs as usize);
        let mut pool = GraphPool::new();
        let gid = pool.alloc(g);
        let cnvtxs = matcher(ctrl, &mut pool, gid);
        assert_eq!(ctrl.mcore, 0, "WCOREPUSH / WCOREPOP unbalanced");
        (pool, gid, cnvtxs)
    }

    #[test]
    fn match_rm_ncon2_island_uses_inclusive_ivecle() {
        // Two islands with vwgt == maxvwgt. ncon > 1 enters the matching
        // branch through ivecle (<=, line 189), so the first island visited
        // takes the next unmatched vertex of perm (lines 192-200).
        let g = graph_from(vec![0, 0, 0], vec![], vec![], 2, vec![2, 2, 2, 2]);
        let mut ctrl = ctrl_with(2, vec![2, 2]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 1);
        assert_eq!(pool.get(gid).cmap, vec![0, 0]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.nvtxs, 1);
        assert_eq!(cg.nedges, 0);
        assert_eq!(cg.vwgt, vec![4, 4]);
    }

    #[test]
    fn match_rm_ncon1_island_uses_strict_bound() {
        // Same shape with ncon == 1: vwgt[i] < maxvwgt[0] is false (line 189),
        // so nothing is matched.
        let g = graph_from(vec![0, 0, 0], vec![], vec![], 1, vec![2, 2]);
        let mut ctrl = ctrl_with(1, vec![2]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 1]);
    }

    #[test]
    fn match_rm_ncon2_pair_blocked_by_second_constraint() {
        // Edge 0-1, vwgt (1,1) each, maxvwgt (2,1): 1*x+y = (2,2) exceeds the
        // second bound (line 225) and 2*x+x = (3,3) fails the 2-hop test
        // (line 232), so each vertex is matched with itself.
        let g = graph_from(vec![0, 1, 2], vec![1, 0], vec![1, 1], 2, vec![1, 1, 1, 1]);
        let mut ctrl = ctrl_with(2, vec![2, 1]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 1]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.nedges, 2);
        assert_eq!(cg.vwgt, vec![1, 1, 1, 1]);
    }

    #[test]
    fn match_rm_ncon2_pair_matches() {
        let g = graph_from(vec![0, 1, 2], vec![1, 0], vec![1, 1], 2, vec![1, 1, 1, 1]);
        let mut ctrl = ctrl_with(2, vec![2, 2]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 1);
        assert_eq!(pool.get(gid).cmap, vec![0, 0]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.nedges, 0);
        assert_eq!(cg.vwgt, vec![2, 2]);
        assert_eq!(cg.tvwgt, vec![2, 2]);
    }

    /// Path 0-1-2 with a heavy centre: the centre cannot match (4 > 3) and
    /// fails the 2-hop test (9 > 3), so it matches itself; each leaf fails
    /// its only edge but passes 2*x+x = 3 <= 3 and is recorded unmatched
    /// (lines 232-235). nunmatched = 2 > 0.1 * 3 triggers Match_2Hop.
    fn heavy_centre_path() -> (graph_t, ctrl_t) {
        let g = graph_from(
            vec![0, 1, 3, 4],
            vec![1, 0, 2, 1],
            vec![1, 1, 1, 1],
            2,
            vec![1, 1, 3, 3, 1, 1],
        );
        (g, ctrl_with(2, vec![3, 3]))
    }

    #[test]
    fn match_rm_ncon2_two_hop_pairs_leaves() {
        let (g, mut ctrl) = heavy_centre_path();
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 1, 0]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.vwgt, vec![2, 2, 3, 3]);
        assert_eq!(cg.nedges, 2);
        assert_eq!(&cg.adjwgt[..2], &[2, 2]);
    }

    #[test]
    fn match_rm_ncon2_no2hop_leaves_unmatched() {
        // Same graph with no2hop: lines 255-264 match the leaves with themselves.
        let (g, mut ctrl) = heavy_centre_path();
        ctrl.no2hop = 1;
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_RM);
        assert_eq!(cnvtxs, 3);
        assert_eq!(pool.get(gid).cmap, vec![0, 1, 2]);
    }

    #[test]
    fn match_shem_ncon2_two_hop_pairs_leaves() {
        // Lines 375-378 and 392-393: the SHEM multi-constraint 2-hop path.
        let (g, mut ctrl) = heavy_centre_path();
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, g, Match_SHEM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 1, 0]);
    }

    /// Path 0-1-2-3 with edge weights 5, 1, 5: every vertex's heaviest
    /// admissible edge belongs to {0,1} or {2,3}, whatever the visit order.
    fn heavy_pairs(ncon: idx_t) -> graph_t {
        graph_from(
            vec![0, 1, 3, 5, 6],
            vec![1, 0, 2, 1, 3, 2],
            vec![5, 5, 1, 1, 5, 5],
            ncon,
            vec![1; (4 * ncon) as usize],
        )
    }

    #[test]
    fn match_shem_ncon2_heavy_edges() {
        let mut ctrl = ctrl_with(2, vec![2, 2]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, heavy_pairs(2), Match_SHEM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 0, 1, 1]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.nedges, 2);
        assert_eq!(&cg.adjncy[..2], &[1, 0]);
        assert_eq!(&cg.adjwgt[..2], &[1, 1]);
        assert_eq!(cg.vwgt, vec![2, 2, 2, 2]);
    }

    #[test]
    fn match_shem_ncon1_heavy_edges() {
        let mut ctrl = ctrl_with(1, vec![2]);
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, heavy_pairs(1), Match_SHEM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 0, 1, 1]);
    }

    #[test]
    fn match_shem_ncon2_tie_uses_better_vbalance() {
        // Centre 0 with adjacency [2, 1], equal edge weights. vwgt: 0 (1,1),
        // 1 (2,0), 2 (1,1); tvwgt (4,2), invtvwgt (0.25,0.5). Visiting 0 first,
        // k = 2 is taken (maxwgt -1 < 1), then k = 1 ties and
        // BetterVBalance(v=(1,1), u1=(1,1), u2=(2,0)) has diff1 = 0.5,
        // diff2 = 0.25, so diff1 - diff2 >= 0 and 1 replaces 2 (lines 365-368).
        let make = || {
            graph_from(
                vec![0, 2, 3, 4],
                vec![2, 1, 0, 0],
                vec![1, 1, 1, 1],
                2,
                vec![1, 1, 2, 0, 1, 1],
            )
        };
        // All three keys are sqrt(1+deg) = 1 (avgdegree 4), so
        // BucketSortKeysInc keeps the irandArrayPermute order: find a seed
        // whose permutation starts at vertex 0.
        let seed = (0u32..1000)
            .find(|&s| {
                let mut st = s;
                let mut tperm = vec![0; 3];
                irandArrayPermute(3, &mut tperm, 0, 1, &mut st);
                tperm[0] == 0
            })
            .expect("a seed visiting vertex 0 first");
        let mut ctrl = ctrl_with(2, vec![3, 3]);
        ctrl.rng_state = seed;
        let (pool, gid, cnvtxs) = run_match(&mut ctrl, make(), Match_SHEM);
        assert_eq!(cnvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 0, 1]);
        let cg = pool.get(pool.get(gid).coarser.expect("coarser linked"));
        assert_eq!(cg.vwgt, vec![3, 1, 1, 1]);
    }

    #[test]
    fn coarsen_graph_ncon2_sets_maxvwgt_and_stops() {
        // tvwgt (2,4), CoarsenTo 1: maxvwgt = 1.5*tvwgt/1 = (3,6) (line 37).
        // After one level nvtxs = 1 is not > CoarsenTo (line 67).
        let g = graph_from(vec![0, 1, 2], vec![1, 0], vec![1, 1], 2, vec![1, 2, 1, 2]);
        let mut ctrl = ctrl_with(2, vec![0, 0]);
        let mut pool = GraphPool::new();
        let gid = pool.alloc(g);
        let cid = CoarsenGraph(&mut ctrl, &mut pool, gid);
        assert_eq!(ctrl.maxvwgt, vec![3, 6]);
        assert_eq!(ctrl.mcore, 0);
        assert_ne!(cid, gid);
        assert_eq!(pool.get(cid).finer, Some(gid));
        assert_eq!(pool.get(cid).nvtxs, 1);
        assert_eq!(pool.get(cid).vwgt, vec![2, 4]);
        assert_eq!(pool.get(gid).cmap, vec![0, 0]);
    }

    #[test]
    fn coarsen_graph_shem_dispatch_on_unequal_weights() {
        // Unequal adjwgt with ctype SHEM: level 1 is Match_SHEM, which pairs
        // {0,1} and {2,3}; the 2-vertex graph passes lines 67-69
        // (2 > 1, 2 < 0.85*4, 2 > 1) and level 2 contracts it to one vertex.
        let mut ctrl = ctrl_with(1, vec![0]);
        ctrl.ctype = METIS_CTYPE_SHEM;
        let mut pool = GraphPool::new();
        let gid = pool.alloc(heavy_pairs(1));
        let cid = CoarsenGraph(&mut ctrl, &mut pool, gid);
        assert_eq!(ctrl.maxvwgt, vec![6]);
        assert_eq!(pool.get(gid).cmap, vec![0, 0, 1, 1]);
        let mid = pool.get(gid).coarser.expect("level 1");
        assert_eq!(pool.get(mid).coarser, Some(cid));
        assert_eq!(pool.get(cid).nvtxs, 1);
        assert_eq!(pool.get(cid).vwgt, vec![4]);
    }

    #[test]
    fn coarsen_graph_nlevels_zero_only_sets_maxvwgt() {
        // tvwgt (5,7), CoarsenTo 3: 7.5/3 = 2.5 -> 2, 10.5/3 = 3.5 -> 3 (line 97).
        let g = graph_from(vec![0, 0, 0, 0], vec![], vec![], 2, vec![2, 1, 2, 3, 1, 3]);
        let mut ctrl = ctrl_with(2, vec![0, 0]);
        ctrl.CoarsenTo = 3;
        let mut pool = GraphPool::new();
        let gid = pool.alloc(g);
        let rid = CoarsenGraphNlevels(&mut ctrl, &mut pool, gid, 0);
        assert_eq!(rid, gid);
        assert_eq!(ctrl.maxvwgt, vec![2, 3]);
        assert!(pool.get(gid).cmap.is_empty());
        assert_eq!(pool.get(gid).coarser, None);
    }

    #[test]
    fn coarsen_graph_nlevels_one_builds_one_level() {
        let mut ctrl = ctrl_with(1, vec![0]);
        ctrl.ctype = METIS_CTYPE_SHEM;
        let mut pool = GraphPool::new();
        let gid = pool.alloc(heavy_pairs(1));
        let rid = CoarsenGraphNlevels(&mut ctrl, &mut pool, gid, 1);
        assert_eq!(pool.get(rid).finer, Some(gid));
        assert_eq!(pool.get(rid).coarser, None);
        assert_eq!(pool.get(rid).nvtxs, 2);
        assert_eq!(pool.get(gid).cmap, vec![0, 0, 1, 1]);
    }
}
