//! HiGHS source: extern/metis/libmetis/balance.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Edge-cut 2-way balancing (used by the edge-based initial separator).
//! `McGeneral2WayBalance` is the `ncon > 1` path, unreachable from HiPO.

use super::{ctrl_t, graph_t, idx_t, real_t, rpq_t};

/// HiGHS: balance.c::Balance2Way (lines 16-34)
pub fn Balance2Way(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t]) {
    // if (ComputeLoadImbalanceDiff(graph, 2, ctrl->pijbm, ctrl->ubfactors) <= 0) return; (lines 18-19)
    if super::mcutil::ComputeLoadImbalanceDiff(graph, 2, &ctrl.pijbm, &ctrl.ubfactors) <= 0.0 {
        return;
    }

    // if (graph->ncon == 1) { (line 21)
    if graph.ncon == 1 {
        // return right away if the balance is OK (lines 22-24)
        // rabs(ntpwgts[0]*graph->tvwgt[0]-graph->pwgts[0]) < 3*graph->tvwgt[0]/graph->nvtxs
        // ntpwgts[0] is real_t (f32), tvwgt[0] is idx_t (i32). C: f32 * i32 -> f32.
        // f32 - i32 -> f32. rabs is fabsf -> f32.
        // 3*graph->tvwgt[0]: int * int -> int. int / int -> int (truncated).
        // Comparison: f32 < int -> int promoted to f32.
        let lhs: f32 = (ntpwgts[0] * (graph.tvwgt[0] as f32) - (graph.pwgts[0] as f32)).abs();
        let rhs: idx_t = (3 * graph.tvwgt[0]) / graph.nvtxs;
        if lhs < (rhs as f32) {
            return;
        }

        // if (graph->nbnd > 0) Bnd2WayBalance(...); else General2WayBalance(...); (lines 26-29)
        if graph.nbnd > 0 {
            Bnd2WayBalance(ctrl, graph, ntpwgts);
        } else {
            General2WayBalance(ctrl, graph, ntpwgts);
        }
    } else {
        // McGeneral2WayBalance(ctrl, graph, ntpwgts); (line 32)
        McGeneral2WayBalance(ctrl, graph, ntpwgts);
    }
}

/// HiGHS: balance.c::Bnd2WayBalance (lines 41-154)
pub fn Bnd2WayBalance(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t]) {
    // WCOREPUSH; (line 51)
    super::wspace::wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj; (line 54)
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy; (line 56)
    // These are read-only pointers in this function; snapshot to avoid borrow conflicts.
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();

    // moved = iwspacemalloc(ctrl, nvtxs); (line 65)
    let mut moved: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // perm = iwspacemalloc(ctrl, nvtxs); (line 66)
    let mut perm: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);

    // Determine from which domain you will be moving data (lines 68-73)
    // tpwgts[0] = graph->tvwgt[0]*ntpwgts[0];
    // tvwgt[0] is idx_t (i32), ntpwgts[0] is real_t (f32). C: i32 * f32 -> f32.
    // Assigned to idx_t (i32), truncating.
    let mut tpwgts: [idx_t; 2] = [0; 2];
    tpwgts[0] = ((graph.tvwgt[0] as f32) * ntpwgts[0]) as idx_t;
    // tpwgts[1] = graph->tvwgt[0] - tpwgts[0];
    tpwgts[1] = graph.tvwgt[0] - tpwgts[0];
    // mindiff = iabs(tpwgts[0]-pwgts[0]);
    let mindiff: idx_t = (tpwgts[0] - graph.pwgts[0]).abs();
    // from = (pwgts[0] < tpwgts[0] ? 1 : 0);
    let from: idx_t = if graph.pwgts[0] < tpwgts[0] { 1 } else { 0 };
    // to = (from+1)%2;
    let to: idx_t = (from + 1) % 2;

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...); (lines 75-78) — debug print, NDEBUG: no-op

    // queue = rpqCreate(nvtxs); (line 80)
    let mut queue: rpq_t = super::pqueue::rpqCreate(nvtxs as usize);

    // iset(nvtxs, -1, moved); (line 82)
    super::gklib::iset(nvtxs as usize, -1, &mut moved);

    // Insert the boundary nodes of the proper partition whose size is OK in the priority queue (lines 85-92)
    let mut nbnd: idx_t = graph.nbnd;
    // irandArrayPermute(nbnd, perm, nbnd/5, 1, &ctrl->rng_state); (line 87)
    super::gklib::irandArrayPermute(nbnd, &mut perm, nbnd / 5, 1, &mut ctrl.rng_state);
    for ii in 0..nbnd as usize {
        let i: idx_t = perm[ii];
        // if (where[bndind[i]] == from && vwgt[bndind[i]] <= mindiff) (line 90)
        let bnd_i: idx_t = graph.bndind[i as usize];
        if graph.where_[bnd_i as usize] == from && graph.vwgt[bnd_i as usize] <= mindiff {
            // rpqInsert(queue, bndind[i], ed[bndind[i]]-id[bndind[i]]); (line 91)
            let key: real_t = (graph.ed[bnd_i as usize] - graph.id[bnd_i as usize]) as real_t;
            super::pqueue::rpqInsert(&mut queue, bnd_i, key);
        }
    }

    let mut mincut: idx_t = graph.mincut;
    // for (nswaps=0; nswaps<nvtxs; nswaps++) { (line 95)
    for nswaps in 0..nvtxs {
        // if ((higain = rpqGetTop(queue)) == -1) break; (lines 96-97)
        let higain: idx_t = super::pqueue::rpqGetTop(&mut queue);
        if higain == -1 {
            break;
        }

        // if (pwgts[to]+vwgt[higain] > tpwgts[to]) break; (lines 99-100)
        if graph.pwgts[to as usize] + graph.vwgt[higain as usize] > tpwgts[to as usize] {
            break;
        }

        // mincut -= (ed[higain]-id[higain]); (line 102)
        mincut -= graph.ed[higain as usize] - graph.id[higain as usize];
        // INC_DEC(pwgts[to], pwgts[from], vwgt[higain]); (line 103)
        // #define INC_DEC(a, b, val) {a += val; b -= val;}
        graph.pwgts[to as usize] += graph.vwgt[higain as usize];
        graph.pwgts[from as usize] -= graph.vwgt[higain as usize];

        // where[higain] = to; (line 105)
        graph.where_[higain as usize] = to;
        // moved[higain] = nswaps; (line 106)
        moved[higain as usize] = nswaps;

        // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...); (lines 108-109) — debug print, NDEBUG: no-op

        // Update the id[i]/ed[i] values of the affected nodes (lines 111-142)
        // SWAP(id[higain], ed[higain], tmp); (line 114)
        std::mem::swap(
            &mut graph.id[higain as usize],
            &mut graph.ed[higain as usize],
        );
        // if (ed[higain] == 0 && xadj[higain] < xadj[higain+1]) BNDDelete(nbnd, bndind, bndptr, higain); (lines 115-116)
        if graph.ed[higain as usize] == 0 && xadj[higain as usize] < xadj[higain as usize + 1] {
            super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }

        // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 118)
        let j_start = xadj[higain as usize] as usize;
        let j_end = xadj[higain as usize + 1] as usize;
        for j in j_start..j_end {
            // k = adjncy[j]; (line 119)
            let k: idx_t = adjncy[j];

            // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 120)
            let kwgt: idx_t = if to == graph.where_[k as usize] {
                graph.adjwgt[j]
            } else {
                -graph.adjwgt[j]
            };
            // INC_DEC(id[k], ed[k], kwgt); (line 121)
            graph.id[k as usize] += kwgt;
            graph.ed[k as usize] -= kwgt;

            // Update its boundary information and queue position (lines 123-141)
            if graph.bndptr[k as usize] != -1 {
                // If k was a boundary vertex
                if graph.ed[k as usize] == 0 {
                    // Not a boundary vertex any more
                    // BNDDelete(nbnd, bndind, bndptr, k); (line 126)
                    super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                    // if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff) rpqDelete(queue, k); (lines 127-128)
                    if moved[k as usize] == -1
                        && graph.where_[k as usize] == from
                        && graph.vwgt[k as usize] <= mindiff
                    {
                        super::pqueue::rpqDelete(&mut queue, k);
                    }
                } else {
                    // If it has not been moved, update its position in the queue
                    // if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff) rpqUpdate(queue, k, ed[k]-id[k]); (lines 131-132)
                    if moved[k as usize] == -1
                        && graph.where_[k as usize] == from
                        && graph.vwgt[k as usize] <= mindiff
                    {
                        let newkey: real_t =
                            (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                        super::pqueue::rpqUpdate(&mut queue, k, newkey);
                    }
                }
            } else {
                if graph.ed[k as usize] > 0 {
                    // It will now become a boundary vertex
                    // BNDInsert(nbnd, bndind, bndptr, k); (line 137)
                    super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                    // if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff) rpqInsert(queue, k, ed[k]-id[k]); (lines 138-139)
                    if moved[k as usize] == -1
                        && graph.where_[k as usize] == from
                        && graph.vwgt[k as usize] <= mindiff
                    {
                        let key: real_t = (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                        super::pqueue::rpqInsert(&mut queue, k, key);
                    }
                }
            }
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...); (lines 145-146) — debug print, NDEBUG: no-op

    // graph->mincut = mincut; (line 148)
    graph.mincut = mincut;
    // graph->nbnd = nbnd; (line 149)
    graph.nbnd = nbnd;

    // rpqDestroy(queue); (line 151)
    super::pqueue::rpqDestroy(queue);

    // WCOREPOP; (line 153)
    super::wspace::wspacepop(ctrl);
}

/// HiGHS: balance.c::General2WayBalance (lines 165-270)
pub fn General2WayBalance(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t]) {
    // WCOREPUSH; (line 175)
    super::wspace::wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj; (line 178)
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy; (line 180)
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();

    // moved = iwspacemalloc(ctrl, nvtxs); (line 189)
    let mut moved: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // perm = iwspacemalloc(ctrl, nvtxs); (line 190)
    let mut perm: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);

    // Determine from which domain you will be moving data (lines 192-197)
    let mut tpwgts: [idx_t; 2] = [0; 2];
    // tpwgts[0] = graph->tvwgt[0]*ntpwgts[0];
    tpwgts[0] = ((graph.tvwgt[0] as f32) * ntpwgts[0]) as idx_t;
    // tpwgts[1] = graph->tvwgt[0] - tpwgts[0];
    tpwgts[1] = graph.tvwgt[0] - tpwgts[0];
    // mindiff = iabs(tpwgts[0]-pwgts[0]);
    let mindiff: idx_t = (tpwgts[0] - graph.pwgts[0]).abs();
    // from = (pwgts[0] < tpwgts[0] ? 1 : 0);
    let from: idx_t = if graph.pwgts[0] < tpwgts[0] { 1 } else { 0 };
    // to = (from+1)%2;
    let to: idx_t = (from + 1) % 2;

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...); (lines 199-201) — debug print, NDEBUG: no-op

    // queue = rpqCreate(nvtxs); (line 203)
    let mut queue: rpq_t = super::pqueue::rpqCreate(nvtxs as usize);

    // iset(nvtxs, -1, moved); (line 205)
    super::gklib::iset(nvtxs as usize, -1, &mut moved);

    // Insert the nodes of the proper partition whose size is OK in the priority queue (lines 208-214)
    // irandArrayPermute(nvtxs, perm, nvtxs/5, 1, &ctrl->rng_state); (line 209)
    super::gklib::irandArrayPermute(nvtxs, &mut perm, nvtxs / 5, 1, &mut ctrl.rng_state);
    for ii in 0..nvtxs as usize {
        let i: idx_t = perm[ii];
        // if (where[i] == from && vwgt[i] <= mindiff) (line 212)
        if graph.where_[i as usize] == from && graph.vwgt[i as usize] <= mindiff {
            // rpqInsert(queue, i, ed[i]-id[i]); (line 213)
            let key: real_t = (graph.ed[i as usize] - graph.id[i as usize]) as real_t;
            super::pqueue::rpqInsert(&mut queue, i, key);
        }
    }

    let mut mincut: idx_t = graph.mincut;
    let mut nbnd: idx_t = graph.nbnd;
    // for (nswaps=0; nswaps<nvtxs; nswaps++) { (line 218)
    for nswaps in 0..nvtxs {
        // if ((higain = rpqGetTop(queue)) == -1) break; (lines 219-220)
        let higain: idx_t = super::pqueue::rpqGetTop(&mut queue);
        if higain == -1 {
            break;
        }

        // if (pwgts[to]+vwgt[higain] > tpwgts[to]) break; (lines 222-223)
        if graph.pwgts[to as usize] + graph.vwgt[higain as usize] > tpwgts[to as usize] {
            break;
        }

        // mincut -= (ed[higain]-id[higain]); (line 225)
        mincut -= graph.ed[higain as usize] - graph.id[higain as usize];
        // INC_DEC(pwgts[to], pwgts[from], vwgt[higain]); (line 226)
        graph.pwgts[to as usize] += graph.vwgt[higain as usize];
        graph.pwgts[from as usize] -= graph.vwgt[higain as usize];

        // where[higain] = to; (line 228)
        graph.where_[higain as usize] = to;
        // moved[higain] = nswaps; (line 229)
        moved[higain as usize] = nswaps;

        // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...); (lines 231-232) — debug print, NDEBUG: no-op

        // Update the id[i]/ed[i] values of the affected nodes (lines 234-258)
        // SWAP(id[higain], ed[higain], tmp); (line 237)
        std::mem::swap(
            &mut graph.id[higain as usize],
            &mut graph.ed[higain as usize],
        );
        // if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1]) BNDDelete(nbnd, bndind, bndptr, higain); (lines 238-239)
        if graph.ed[higain as usize] == 0
            && graph.bndptr[higain as usize] != -1
            && xadj[higain as usize] < xadj[higain as usize + 1]
        {
            super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }
        // if (ed[higain] > 0 && bndptr[higain] == -1) BNDInsert(nbnd, bndind, bndptr, higain); (lines 240-241)
        if graph.ed[higain as usize] > 0 && graph.bndptr[higain as usize] == -1 {
            super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }

        // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 243)
        let j_start = xadj[higain as usize] as usize;
        let j_end = xadj[higain as usize + 1] as usize;
        for j in j_start..j_end {
            // k = adjncy[j]; (line 244)
            let k: idx_t = adjncy[j];

            // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 246)
            let kwgt: idx_t = if to == graph.where_[k as usize] {
                graph.adjwgt[j]
            } else {
                -graph.adjwgt[j]
            };
            // INC_DEC(id[k], ed[k], kwgt); (line 247)
            graph.id[k as usize] += kwgt;
            graph.ed[k as usize] -= kwgt;

            // Update the queue position (lines 249-251)
            // if (moved[k] == -1 && where[k] == from && vwgt[k] <= mindiff) rpqUpdate(queue, k, ed[k]-id[k]);
            if moved[k as usize] == -1
                && graph.where_[k as usize] == from
                && graph.vwgt[k as usize] <= mindiff
            {
                let newkey: real_t = (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                super::pqueue::rpqUpdate(&mut queue, k, newkey);
            }

            // Update its boundary information (lines 253-257)
            // if (ed[k] == 0 && bndptr[k] != -1) BNDDelete(nbnd, bndind, bndptr, k);
            if graph.ed[k as usize] == 0 && graph.bndptr[k as usize] != -1 {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
            } else if graph.ed[k as usize] > 0 && graph.bndptr[k as usize] == -1 {
                // else if (ed[k] > 0 && bndptr[k] == -1) BNDInsert(nbnd, bndind, bndptr, k);
                super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
            }
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...); (lines 261-262) — debug print, NDEBUG: no-op

    // graph->mincut = mincut; (line 264)
    graph.mincut = mincut;
    // graph->nbnd = nbnd; (line 265)
    graph.nbnd = nbnd;

    // rpqDestroy(queue); (line 267)
    super::pqueue::rpqDestroy(queue);

    // WCOREPOP; (line 269)
    super::wspace::wspacepop(ctrl);
}

/// HiGHS: balance.c::McGeneral2WayBalance (lines 276-490)
/// `ncon > 1` only; never reached from HiPO (`Highs_METIS_NodeND` uses ncon 1).
pub fn McGeneral2WayBalance(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t]) {
    // WCOREPUSH; (line 287) — macros.h: `wspacepush(ctrl)`
    super::wspace::wspacepush(ctrl);

    let nvtxs = graph.nvtxs;
    let ncon = graph.ncon;
    // xadj / adjncy: `graph->xadj_orig ? graph->xadj_orig : graph->xadj` (never written here)
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj; (line 291)
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy; (line 293)
    // vwgt, adjwgt, invtvwgt, where, id, ed, pwgts, bndptr, bndind from graph (lines 292-301)

    // moved = iwspacemalloc(ctrl, nvtxs); (line 303)
    let mut moved = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // swaps = iwspacemalloc(ctrl, nvtxs); (line 304)
    let mut swaps = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // perm = iwspacemalloc(ctrl, nvtxs); (line 305)
    let mut perm = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // qnum = iwspacemalloc(ctrl, nvtxs); (line 306)
    let mut qnum = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // newbalv = rwspacemalloc(ctrl, ncon); (line 307)
    let mut newbalv = super::wspace::rwspacemalloc(ctrl, ncon);
    // minbalv = rwspacemalloc(ctrl, ncon); (line 308)
    let mut minbalv = super::wspace::rwspacemalloc(ctrl, ncon);
    // qsizes = iwspacemalloc(ctrl, 2*ncon); (line 309)
    let mut qsizes = super::wspace::iwspacemalloc(ctrl, 2 * ncon);

    // limit = gk_min(gk_max(0.01*nvtxs, 15), 100); (line 311)
    // gk_max(a, b) is `((a) >= (b) ? (a) : (b))`, gk_min(a, b) is `((a) >= (b) ? (b) : (a))`.
    // 0.01*nvtxs: 0.01 is a double literal, nvtxs is idx_t (i32). C promotes nvtxs to double,
    // multiplies in double, then gk_max compares double vs int 15 (int promoted to double).
    // Result is double. gk_min compares that double vs int 100 (promoted to double).
    // limit is idx_t, so the final result is truncated to i32.
    let tmp_max: f64 = if (0.01 * (nvtxs as f64)) >= 15.0 {
        0.01 * (nvtxs as f64)
    } else {
        15.0
    };
    let limit: idx_t = if tmp_max >= 100.0 {
        100
    } else {
        tmp_max as idx_t
    };

    // queues = (rpq_t **)wspacemalloc(ctrl, 2*ncon*sizeof(rpq_t *)); (line 314)
    // for (i=0; i<2*ncon; i++) { queues[i] = rpqCreate(nvtxs); qsizes[i] = 0; } (lines 315-318)
    let mut queues: Vec<rpq_t> = Vec::with_capacity((2 * ncon) as usize);
    for i in 0..(2 * ncon) as usize {
        queues.push(super::pqueue::rpqCreate(nvtxs as usize));
        qsizes[i] = 0;
    }

    // for (i=0; i<nvtxs; i++) { qnum[i] = iargmax_nrm(ncon, vwgt+i*ncon, invtvwgt); qsizes[2*qnum[i]+where[i]]++; } (lines 320-323)
    for i in 0..nvtxs as usize {
        qnum[i] = super::mcutil::iargmax_nrm(
            ncon as usize,
            &graph.vwgt[i * ncon as usize..(i + 1) * ncon as usize],
            &graph.invtvwgt,
        );
        qsizes[(2 * qnum[i] + graph.where_[i]) as usize] += 1;
    }

    // for the empty queues, move into them vertices from other queues (lines 327-345)
    for from in 0..2 {
        for j in 0..ncon {
            if qsizes[(2 * j + from) as usize] == 0 {
                for i in 0..nvtxs as usize {
                    if graph.where_[i] != from {
                        continue;
                    }

                    let k = super::mcutil::iargmax2_nrm(
                        ncon as usize,
                        &graph.vwgt[i * ncon as usize..(i + 1) * ncon as usize],
                        &graph.invtvwgt,
                    );
                    // if (k == j &&
                    //     qsizes[2*qnum[i]+from] > qsizes[2*j+from] &&
                    //     vwgt[i*ncon+qnum[i]]*invtvwgt[qnum[i]] < 1.3*vwgt[i*ncon+j]*invtvwgt[j])
                    // (lines 335-337)
                    // vwgt is idx_t (i32), invtvwgt is real_t (f32).
                    // vwgt[i*ncon+qnum[i]]*invtvwgt[qnum[i]]: i32 * f32 -> f32 (C: int promoted to float)
                    // 1.3*vwgt[i*ncon+j]*invtvwgt[j]: 1.3 is double, vwgt is i32 -> double * i32 -> double,
                    //   then double * f32 -> double.
                    // Comparison: f32 < double -> f32 promoted to double.
                    if k == j
                        && qsizes[(2 * qnum[i] + from) as usize] > qsizes[(2 * j + from) as usize]
                        && (((graph.vwgt[i * ncon as usize + qnum[i] as usize] as f32)
                            * graph.invtvwgt[qnum[i] as usize]) as f64)
                            < (1.3 * (graph.vwgt[i * ncon as usize + j as usize] as f64))
                                * (graph.invtvwgt[j as usize] as f64)
                    {
                        qsizes[(2 * qnum[i] + from) as usize] -= 1;
                        qsizes[(2 * j + from) as usize] += 1;
                        qnum[i] = j;
                    }
                }
            }
        }
    }

    // minbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl->pijbm, ctrl->ubfactors, minbalv); (line 348)
    let mut minbal = super::mcutil::ComputeLoadImbalanceDiffVec(
        graph,
        2,
        &ctrl.pijbm,
        &ctrl.ubfactors,
        &mut minbalv,
    );

    let mut newcut = graph.mincut;
    let mut mincut = graph.mincut;
    let mut mincutorder: idx_t = -1;

    // if (ctrl->dbglvl&METIS_DBG_REFINE) { ... } (lines 353-360) — debug print, NDEBUG: no-op

    // iset(nvtxs, -1, moved); (line 362)
    super::gklib::iset(nvtxs as usize, -1, &mut moved);

    // nbnd = graph->nbnd; (line 366)
    let mut nbnd = graph.nbnd;

    // irandArrayPermute(nvtxs, perm, nvtxs/10, 1, &ctrl->rng_state); (line 367)
    super::gklib::irandArrayPermute(nvtxs, &mut perm, nvtxs / 10, 1, &mut ctrl.rng_state);

    // for (ii=0; ii<nvtxs; ii++) { i = perm[ii]; rpqInsert(queues[2*qnum[i]+where[i]], i, ed[i]-id[i]); } (lines 368-371)
    for ii in 0..nvtxs as usize {
        let i = perm[ii];
        let qi = (2 * qnum[i as usize] + graph.where_[i as usize]) as usize;
        let key = graph.ed[i as usize] - graph.id[i as usize];
        super::pqueue::rpqInsert(&mut queues[qi], i, key as real_t);
    }

    // for (nswaps=0; nswaps<nvtxs; nswaps++) { (line 373)
    let mut nswaps: idx_t = 0;
    while nswaps < nvtxs {
        // if (minbal <= 0.0) break; (lines 374-375)
        if minbal <= 0.0 {
            break;
        }

        // SelectQueue(graph, ctrl->pijbm, ctrl->ubfactors, queues, &from, &cnum); (line 377)
        let mut from: idx_t = 0;
        let mut cnum: idx_t = 0;
        super::fm::SelectQueue(
            graph,
            &ctrl.pijbm,
            &ctrl.ubfactors,
            &queues,
            &mut from,
            &mut cnum,
        );

        // to = (from+1)%2; (line 378)
        let to = (from + 1) % 2;

        // if (from == -1 || (higain = rpqGetTop(queues[2*cnum+from])) == -1) break; (lines 380-381)
        if from == -1 {
            break;
        }
        let higain = super::pqueue::rpqGetTop(&mut queues[(2 * cnum + from) as usize]);
        if higain == -1 {
            break;
        }

        // newcut -= (ed[higain]-id[higain]); (line 383)
        newcut -= graph.ed[higain as usize] - graph.id[higain as usize];

        // iaxpy(ncon,  1, vwgt+higain*ncon, 1, pwgts+to*ncon,   1); (line 385)
        super::gklib::iaxpy(
            ncon as usize,
            1,
            &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
            1,
            &mut graph.pwgts[to as usize * ncon as usize..(to as usize + 1) * ncon as usize],
            1,
        );
        // iaxpy(ncon, -1, vwgt+higain*ncon, 1, pwgts+from*ncon, 1); (line 386)
        super::gklib::iaxpy(
            ncon as usize,
            -1,
            &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
            1,
            &mut graph.pwgts[from as usize * ncon as usize..(from as usize + 1) * ncon as usize],
            1,
        );

        // newbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl->pijbm, ctrl->ubfactors, newbalv); (line 387)
        let newbal = super::mcutil::ComputeLoadImbalanceDiffVec(
            graph,
            2,
            &ctrl.pijbm,
            &ctrl.ubfactors,
            &mut newbalv,
        );

        // if (newbal < minbal || (newbal == minbal &&
        //     (newcut < mincut ||
        //      (newcut == mincut && BetterBalance2Way(ncon, minbalv, newbalv))))) (lines 389-391)
        if newbal < minbal
            || (newbal == minbal
                && (newcut < mincut
                    || (newcut == mincut
                        && super::mcutil::BetterBalance2Way(ncon, &minbalv, &newbalv) != 0)))
        {
            mincut = newcut;
            minbal = newbal;
            mincutorder = nswaps;
            super::gklib::rcopy(ncon as usize, &newbalv, &mut minbalv);
        } else if nswaps - mincutorder > limit {
            // We hit the limit, undo last move (lines 397-402)
            // C: newcut += (ed[higain]-id[higain]); newcut is dead after the loop, so the store is elided.
            super::gklib::iaxpy(
                ncon as usize,
                1,
                &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
                1,
                &mut graph.pwgts
                    [from as usize * ncon as usize..(from as usize + 1) * ncon as usize],
                1,
            );
            super::gklib::iaxpy(
                ncon as usize,
                -1,
                &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
                1,
                &mut graph.pwgts[to as usize * ncon as usize..(to as usize + 1) * ncon as usize],
                1,
            );
            break;
        }

        // where[higain] = to; (line 404)
        graph.where_[higain as usize] = to;
        // moved[higain] = nswaps; (line 405)
        moved[higain as usize] = nswaps;
        // swaps[nswaps] = higain; (line 406)
        swaps[nswaps as usize] = higain;

        // if (ctrl->dbglvl&METIS_DBG_MOVEINFO) { ... } (lines 408-414) — debug print, NDEBUG: no-op

        // SWAP(id[higain], ed[higain], tmp); (line 420)
        std::mem::swap(
            &mut graph.id[higain as usize],
            &mut graph.ed[higain as usize],
        );

        // if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        //   BNDDelete(nbnd, bndind, bndptr, higain); (lines 421-422)
        if graph.ed[higain as usize] == 0
            && graph.bndptr[higain as usize] != -1
            && xadj[higain as usize] < xadj[higain as usize + 1]
        {
            super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }
        // if (ed[higain] > 0 && bndptr[higain] == -1)
        //   BNDInsert(nbnd, bndind, bndptr, higain); (lines 423-424)
        if graph.ed[higain as usize] > 0 && graph.bndptr[higain as usize] == -1 {
            super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }

        // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 426)
        let j_start = xadj[higain as usize] as usize;
        let j_end = xadj[higain as usize + 1] as usize;
        for j in j_start..j_end {
            // k = adjncy[j]; (line 427)
            let k = adjncy[j] as usize;

            // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 429)
            let kwgt = if to == graph.where_[k] {
                graph.adjwgt[j]
            } else {
                -graph.adjwgt[j]
            };

            // INC_DEC(id[k], ed[k], kwgt); (line 430)
            // #define INC_DEC(a, b, val) {a += val; b -= val;}
            graph.id[k] += kwgt;
            graph.ed[k] -= kwgt;

            // if (moved[k] == -1)
            //   rpqUpdate(queues[2*qnum[k]+where[k]], k, ed[k]-id[k]); (lines 433-434)
            if moved[k] == -1 {
                let qi = (2 * qnum[k] + graph.where_[k]) as usize;
                let newkey = (graph.ed[k] - graph.id[k]) as real_t;
                super::pqueue::rpqUpdate(&mut queues[qi], k as idx_t, newkey);
            }

            // if (ed[k] == 0 && bndptr[k] != -1)
            //   BNDDelete(nbnd, bndind, bndptr, k); (lines 437-438)
            if graph.ed[k] == 0 && graph.bndptr[k] != -1 {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k as idx_t);
            } else if graph.ed[k] > 0 && graph.bndptr[k] == -1 {
                // else if (ed[k] > 0 && bndptr[k] == -1)
                //   BNDInsert(nbnd, bndind, bndptr, k); (lines 439-440)
                super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k as idx_t);
            }
        }

        nswaps += 1;
    }

    // Roll back computations (lines 446-472)
    // for (nswaps--; nswaps>mincutorder; nswaps--) { (line 449)
    nswaps -= 1;
    while nswaps > mincutorder {
        // higain = swaps[nswaps]; (line 450)
        let higain = swaps[nswaps as usize];

        // to = where[higain] = (where[higain]+1)%2; (line 452)
        let to = (graph.where_[higain as usize] + 1) % 2;
        graph.where_[higain as usize] = to;

        // SWAP(id[higain], ed[higain], tmp); (line 453)
        std::mem::swap(
            &mut graph.id[higain as usize],
            &mut graph.ed[higain as usize],
        );

        // if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
        //   BNDDelete(nbnd, bndind, bndptr, higain); (lines 454-455)
        if graph.ed[higain as usize] == 0
            && graph.bndptr[higain as usize] != -1
            && xadj[higain as usize] < xadj[higain as usize + 1]
        {
            super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        } else if graph.ed[higain as usize] > 0 && graph.bndptr[higain as usize] == -1 {
            // else if (ed[higain] > 0 && bndptr[higain] == -1)
            //   BNDInsert(nbnd, bndind, bndptr, higain); (lines 456-457)
            super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        }

        // iaxpy(ncon,  1, vwgt+higain*ncon, 1, pwgts+to*ncon,         1); (line 459)
        super::gklib::iaxpy(
            ncon as usize,
            1,
            &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
            1,
            &mut graph.pwgts[to as usize * ncon as usize..(to as usize + 1) * ncon as usize],
            1,
        );
        // iaxpy(ncon, -1, vwgt+higain*ncon, 1, pwgts+((to+1)%2)*ncon, 1); (line 460)
        let other = (to + 1) % 2;
        super::gklib::iaxpy(
            ncon as usize,
            -1,
            &graph.vwgt[higain as usize * ncon as usize..(higain as usize + 1) * ncon as usize],
            1,
            &mut graph.pwgts[other as usize * ncon as usize..(other as usize + 1) * ncon as usize],
            1,
        );

        // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 461)
        let j_start = xadj[higain as usize] as usize;
        let j_end = xadj[higain as usize + 1] as usize;
        for j in j_start..j_end {
            // k = adjncy[j]; (line 462)
            let k = adjncy[j] as usize;

            // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 464)
            let kwgt = if to == graph.where_[k] {
                graph.adjwgt[j]
            } else {
                -graph.adjwgt[j]
            };

            // INC_DEC(id[k], ed[k], kwgt); (line 465)
            graph.id[k] += kwgt;
            graph.ed[k] -= kwgt;

            // if (bndptr[k] != -1 && ed[k] == 0)
            //   BNDDelete(nbnd, bndind, bndptr, k); (lines 467-468)
            if graph.bndptr[k] != -1 && graph.ed[k] == 0 {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k as idx_t);
            }
            // if (bndptr[k] == -1 && ed[k] > 0)
            //   BNDInsert(nbnd, bndind, bndptr, k); (lines 469-470)
            if graph.bndptr[k] == -1 && graph.ed[k] > 0 {
                super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k as idx_t);
            }
        }

        nswaps -= 1;
    }

    // if (ctrl->dbglvl&METIS_DBG_REFINE) { ... } (lines 474-480) — debug print, NDEBUG: no-op

    // graph->mincut = mincut; (line 482)
    graph.mincut = mincut;
    // graph->nbnd = nbnd; (line 483)
    graph.nbnd = nbnd;

    // for (i=0; i<2*ncon; i++) rpqDestroy(queues[i]); (lines 486-487)
    for queue in queues {
        super::pqueue::rpqDestroy(queue);
    }

    // WCOREPOP; (line 489) — macros.h: `wspacepop(ctrl)`
    super::wspace::wspacepop(ctrl);
}
