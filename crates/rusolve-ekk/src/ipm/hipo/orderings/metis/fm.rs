//! HiGHS source: extern/metis/libmetis/fm.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Edge-cut FM refinement (`FM_2WayCutRefine`) for the edge-based initial
//! separator; the `Mc` variant and `SelectQueue` serve `ncon > 1` only.

use super::{ctrl_t, graph_t, idx_t, real_t, rpq_t};

/// HiGHS: fm.c::FM_2WayRefine (lines 17-23)
pub fn FM_2WayRefine(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t], niter: idx_t) {
    // if (graph->ncon == 1) FM_2WayCutRefine(...); else FM_Mc2WayCutRefine(...); (lines 19-22)
    if graph.ncon == 1 {
        FM_2WayCutRefine(ctrl, graph, ntpwgts, niter);
    } else {
        FM_Mc2WayCutRefine(ctrl, graph, ntpwgts, niter);
    }
}

/// HiGHS: fm.c::FM_2WayCutRefine (lines 29-197)
pub fn FM_2WayCutRefine(ctrl: &mut ctrl_t, graph: &mut graph_t, ntpwgts: &[real_t], niter: idx_t) {
    // WCOREPUSH; (line 39)
    super::wspace::wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    // xadj = graph->xadj_orig ? graph->xadj_orig : graph->xadj; (line 42)
    // adjncy = graph->adjncy_orig ? graph->adjncy_orig : graph->adjncy; (line 44)
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();

    // moved = iwspacemalloc(ctrl, nvtxs); (line 53)
    let mut moved: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // swaps = iwspacemalloc(ctrl, nvtxs); (line 54)
    let mut swaps: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    // perm = iwspacemalloc(ctrl, nvtxs); (line 55)
    let mut perm: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);

    // tpwgts[0] = graph->tvwgt[0]*ntpwgts[0]; (line 57)
    // tvwgt[0] is idx_t (i32), ntpwgts[0] is real_t (f32). C: i32 * f32 -> f32.
    // Assigned to idx_t (i32), truncating.
    let mut tpwgts: [idx_t; 2] = [0; 2];
    tpwgts[0] = ((graph.tvwgt[0] as f32) * ntpwgts[0]) as idx_t;
    // tpwgts[1] = graph->tvwgt[0]-tpwgts[0]; (line 58)
    tpwgts[1] = graph.tvwgt[0] - tpwgts[0];

    // limit = gk_min(gk_max(0.01*nvtxs, 15), 100); (line 60)
    // gk_max(a,b) is (a >= b ? a : b). gk_min(a,b) is (a >= b ? b : a).
    // 0.01*nvtxs is f64. 15 is int, promoted to f64. Result f64.
    // gk_max(0.01*nvtxs, 15.0) -> f64.
    // gk_min(that, 100.0) -> f64.
    // Assigned to idx_t (i32), truncating.
    let tmp_max: f64 = if (0.01 * (nvtxs as f64)) >= 15.0 {
        0.01 * (nvtxs as f64)
    } else {
        15.0
    };
    let limit: idx_t = (if tmp_max >= 100.0 { 100.0 } else { tmp_max }) as idx_t;

    // avgvwgt = gk_min((pwgts[0]+pwgts[1])/20, 2*(pwgts[0]+pwgts[1])/nvtxs); (line 61)
    // All operands are idx_t (i32). Integer division.
    // gk_min(a,b) is (a >= b ? b : a).
    let a1: idx_t = (graph.pwgts[0] + graph.pwgts[1]) / 20;
    let a2: idx_t = 2 * (graph.pwgts[0] + graph.pwgts[1]) / nvtxs;
    let avgvwgt: idx_t = if a1 >= a2 { a2 } else { a1 };

    // queues[0] = rpqCreate(nvtxs); (line 63)
    // queues[1] = rpqCreate(nvtxs); (line 64)
    let mut queues: [rpq_t; 2] = [
        super::pqueue::rpqCreate(nvtxs as usize),
        super::pqueue::rpqCreate(nvtxs as usize),
    ];

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, Print2WayRefineStats(...)); (lines 66-67)
    // Dropped: debug-only under NDEBUG.

    // origdiff = iabs(tpwgts[0]-pwgts[0]); (line 69)
    let origdiff: idx_t = (tpwgts[0] - graph.pwgts[0]).abs();
    // iset(nvtxs, -1, moved); (line 70)
    super::gklib::iset(nvtxs as usize, -1, &mut moved);

    // for (pass=0; pass<niter; pass++) { (line 71)
    let mut pass: idx_t = 0;
    while pass < niter {
        // rpqReset(queues[0]); (line 72)
        super::pqueue::rpqReset(&mut queues[0]);
        // rpqReset(queues[1]); (line 73)
        super::pqueue::rpqReset(&mut queues[1]);

        let mut mincutorder: idx_t = -1;
        // newcut = mincut = initcut = graph->mincut; (line 76)
        let mut newcut: idx_t = graph.mincut;
        let mut mincut: idx_t = graph.mincut;
        let initcut: idx_t = graph.mincut;
        // mindiff = iabs(tpwgts[0]-pwgts[0]); (line 77)
        let mut mindiff: idx_t = (tpwgts[0] - graph.pwgts[0]).abs();

        // Insert boundary nodes in the priority queues (lines 80-86)
        let mut nbnd: idx_t = graph.nbnd;
        // irandArrayPermute(nbnd, perm, nbnd, 1, &ctrl->rng_state); (line 82)
        super::gklib::irandArrayPermute(nbnd, &mut perm, nbnd, 1, &mut ctrl.rng_state);
        for ii in 0..nbnd as usize {
            let i: idx_t = perm[ii];
            // rpqInsert(queues[where[bndind[i]]], bndind[i], ed[bndind[i]]-id[bndind[i]]); (line 85)
            let bnd_i: idx_t = graph.bndind[i as usize];
            let qidx: usize = graph.where_[bnd_i as usize] as usize;
            let key: real_t = (graph.ed[bnd_i as usize] - graph.id[bnd_i as usize]) as real_t;
            super::pqueue::rpqInsert(&mut queues[qidx], bnd_i, key);
        }

        // for (nswaps=0; nswaps<nvtxs; nswaps++) { (line 88)
        let mut nswaps: idx_t = 0;
        while nswaps < nvtxs {
            // from = (tpwgts[0]-pwgts[0] < tpwgts[1]-pwgts[1] ? 0 : 1); (line 89)
            let from: idx_t = if (tpwgts[0] - graph.pwgts[0]) < (tpwgts[1] - graph.pwgts[1]) {
                0
            } else {
                1
            };
            // to = (from+1)%2; (line 90)
            let to: idx_t = (from + 1) % 2;

            // if ((higain = rpqGetTop(queues[from])) == -1) break; (lines 92-93)
            let higain: idx_t = super::pqueue::rpqGetTop(&mut queues[from as usize]);
            if higain == -1 {
                break;
            }

            // newcut -= (ed[higain]-id[higain]); (line 95)
            newcut -= graph.ed[higain as usize] - graph.id[higain as usize];
            // INC_DEC(pwgts[to], pwgts[from], vwgt[higain]); (line 96)
            graph.pwgts[to as usize] += graph.vwgt[higain as usize];
            graph.pwgts[from as usize] -= graph.vwgt[higain as usize];

            // if ((newcut < mincut && iabs(tpwgts[0]-pwgts[0]) <= origdiff+avgvwgt) ||
            //     (newcut == mincut && iabs(tpwgts[0]-pwgts[0]) < mindiff)) { (lines 98-99)
            let abs_diff: idx_t = (tpwgts[0] - graph.pwgts[0]).abs();
            if (newcut < mincut && abs_diff <= origdiff + avgvwgt)
                || (newcut == mincut && abs_diff < mindiff)
            {
                // mincut = newcut; (line 100)
                mincut = newcut;
                // mindiff = iabs(tpwgts[0]-pwgts[0]); (line 101)
                mindiff = abs_diff;
                // mincutorder = nswaps; (line 102)
                mincutorder = nswaps;
            } else if nswaps - mincutorder > limit {
                // We hit the limit, undo last move (lines 104-108)
                // newcut += (ed[higain]-id[higain]); (line 105)
                newcut += graph.ed[higain as usize] - graph.id[higain as usize];
                // PORT NOTE: C leaves this store dead (the loop breaks and newcut is reset each pass);
                // read it once so the Rust unused_assignments lint stays quiet. No semantic effect.
                let _ = newcut;
                // INC_DEC(pwgts[from], pwgts[to], vwgt[higain]); (line 106)
                graph.pwgts[from as usize] += graph.vwgt[higain as usize];
                graph.pwgts[to as usize] -= graph.vwgt[higain as usize];
                // break; (line 107)
                break;
            }

            // where[higain] = to; (line 110)
            graph.where_[higain as usize] = to;
            // moved[higain] = nswaps; (line 111)
            moved[higain as usize] = nswaps;
            // swaps[nswaps] = higain; (line 112)
            swaps[nswaps as usize] = higain;

            // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...); (lines 114-115)
            // Dropped: debug-only under NDEBUG.

            // Update the id[i]/ed[i] values of the affected nodes (lines 117-149)
            // SWAP(id[higain], ed[higain], tmp); (line 120)
            std::mem::swap(
                &mut graph.id[higain as usize],
                &mut graph.ed[higain as usize],
            );
            // if (ed[higain] == 0 && xadj[higain] < xadj[higain+1]) BNDDelete(nbnd, bndind, bndptr, higain); (lines 121-122)
            if graph.ed[higain as usize] == 0 && xadj[higain as usize] < xadj[higain as usize + 1] {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            }

            // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 124)
            let j_start = xadj[higain as usize] as usize;
            let j_end = xadj[higain as usize + 1] as usize;
            for j in j_start..j_end {
                // k = adjncy[j]; (line 125)
                let k: idx_t = adjncy[j];

                // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 127)
                let kwgt: idx_t = if to == graph.where_[k as usize] {
                    graph.adjwgt[j]
                } else {
                    -graph.adjwgt[j]
                };
                // INC_DEC(id[k], ed[k], kwgt); (line 128)
                graph.id[k as usize] += kwgt;
                graph.ed[k as usize] -= kwgt;

                // Update its boundary information and queue position (lines 130-148)
                if graph.bndptr[k as usize] != -1 {
                    // If k was a boundary vertex
                    if graph.ed[k as usize] == 0 {
                        // Not a boundary vertex any more
                        // BNDDelete(nbnd, bndind, bndptr, k); (line 133)
                        super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                        // if (moved[k] == -1) rpqDelete(queues[where[k]], k); (lines 134-135)
                        if moved[k as usize] == -1 {
                            let qidx: usize = graph.where_[k as usize] as usize;
                            super::pqueue::rpqDelete(&mut queues[qidx], k);
                        }
                    } else {
                        // If it has not been moved, update its position in the queue
                        // if (moved[k] == -1) rpqUpdate(queues[where[k]], k, ed[k]-id[k]); (lines 138-139)
                        if moved[k as usize] == -1 {
                            let qidx: usize = graph.where_[k as usize] as usize;
                            let newkey: real_t =
                                (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                            super::pqueue::rpqUpdate(&mut queues[qidx], k, newkey);
                        }
                    }
                } else {
                    if graph.ed[k as usize] > 0 {
                        // It will now become a boundary vertex
                        // BNDInsert(nbnd, bndind, bndptr, k); (line 144)
                        super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                        // if (moved[k] == -1) rpqInsert(queues[where[k]], k, ed[k]-id[k]); (lines 145-146)
                        if moved[k as usize] == -1 {
                            let qidx: usize = graph.where_[k as usize] as usize;
                            let key: real_t =
                                (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                            super::pqueue::rpqInsert(&mut queues[qidx], k, key);
                        }
                    }
                }
            }

            nswaps += 1;
        }

        // Roll back computations (lines 154-181)
        // for (i=0; i<nswaps; i++) moved[swaps[i]] = -1; (lines 157-158)
        for i in 0..nswaps as usize {
            moved[swaps[i] as usize] = -1;
        }
        // for (nswaps--; nswaps>mincutorder; nswaps--) { (line 159)
        nswaps -= 1;
        while nswaps > mincutorder {
            // higain = swaps[nswaps]; (line 160)
            let higain: idx_t = swaps[nswaps as usize];

            // to = where[higain] = (where[higain]+1)%2; (line 162)
            let to: idx_t = (graph.where_[higain as usize] + 1) % 2;
            graph.where_[higain as usize] = to;

            // SWAP(id[higain], ed[higain], tmp); (line 163)
            std::mem::swap(
                &mut graph.id[higain as usize],
                &mut graph.ed[higain as usize],
            );
            // if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1]) BNDDelete(nbnd, bndind, bndptr, higain); (lines 164-165)
            if graph.ed[higain as usize] == 0
                && graph.bndptr[higain as usize] != -1
                && xadj[higain as usize] < xadj[higain as usize + 1]
            {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            } else if graph.ed[higain as usize] > 0 && graph.bndptr[higain as usize] == -1 {
                // else if (ed[higain] > 0 && bndptr[higain] == -1) BNDInsert(nbnd, bndind, bndptr, higain); (lines 166-167)
                super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            }

            // INC_DEC(pwgts[to], pwgts[(to+1)%2], vwgt[higain]); (line 169)
            let other: idx_t = (to + 1) % 2;
            graph.pwgts[to as usize] += graph.vwgt[higain as usize];
            graph.pwgts[other as usize] -= graph.vwgt[higain as usize];

            // for (j=xadj[higain]; j<xadj[higain+1]; j++) { (line 170)
            let j_start = xadj[higain as usize] as usize;
            let j_end = xadj[higain as usize + 1] as usize;
            for j in j_start..j_end {
                // k = adjncy[j]; (line 171)
                let k: idx_t = adjncy[j];

                // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 173)
                let kwgt: idx_t = if to == graph.where_[k as usize] {
                    graph.adjwgt[j]
                } else {
                    -graph.adjwgt[j]
                };
                // INC_DEC(id[k], ed[k], kwgt); (line 174)
                graph.id[k as usize] += kwgt;
                graph.ed[k as usize] -= kwgt;

                // if (bndptr[k] != -1 && ed[k] == 0) BNDDelete(nbnd, bndind, bndptr, k); (lines 176-177)
                if graph.bndptr[k as usize] != -1 && graph.ed[k as usize] == 0 {
                    super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                }
                // if (bndptr[k] == -1 && ed[k] > 0) BNDInsert(nbnd, bndind, bndptr, k); (lines 178-179)
                if graph.bndptr[k as usize] == -1 && graph.ed[k as usize] > 0 {
                    super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                }
            }

            nswaps -= 1;
        }

        // graph->mincut = mincut; (line 183)
        graph.mincut = mincut;
        // graph->nbnd = nbnd; (line 184)
        graph.nbnd = nbnd;

        // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, Print2WayRefineStats(...)); (lines 186-187)
        // Dropped: debug-only under NDEBUG.

        // if (mincutorder <= 0 || mincut == initcut) break; (lines 189-190)
        if mincutorder <= 0 || mincut == initcut {
            break;
        }

        pass += 1;
    }

    // rpqDestroy(queues[0]); (line 193)
    let [q0, q1] = queues;
    super::pqueue::rpqDestroy(q0);
    // rpqDestroy(queues[1]); (line 194)
    super::pqueue::rpqDestroy(q1);

    // WCOREPOP; (line 196)
    super::wspace::wspacepop(ctrl);
}

/// HiGHS: fm.c::FM_Mc2WayCutRefine (lines 203-424)
/// `ncon > 1` only; never reached from HiPO.
pub fn FM_Mc2WayCutRefine(
    ctrl: &mut ctrl_t,
    graph: &mut graph_t,
    ntpwgts: &[real_t],
    niter: idx_t,
) {
    // PORT NOTE: ntpwgts is unused in the body (only passed to the dropped debug printer).
    let _ = ntpwgts;

    // WCOREPUSH (line 216)
    super::wspace::wspacepush(ctrl);

    let nvtxs: idx_t = graph.nvtxs;
    let ncon: idx_t = graph.ncon;

    // PORT NOTE: xadj and adjncy are C pointers to xadj_orig/adjncy_orig when set.
    // The rest of the function mutates graph fields but never the xadj/adjncy contents,
    // so we take an immutable snapshot to avoid borrow conflicts.
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();

    // Workspace allocations (lines 232-238)
    let mut moved: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    let mut swaps: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    let mut perm: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    let mut qnum: Vec<idx_t> = super::wspace::iwspacemalloc(ctrl, nvtxs);
    let mut ubfactors: Vec<real_t> = super::wspace::rwspacemalloc(ctrl, ncon);
    let mut newbalv: Vec<real_t> = super::wspace::rwspacemalloc(ctrl, ncon);
    let mut minbalv: Vec<real_t> = super::wspace::rwspacemalloc(ctrl, ncon);

    // limit = gk_min(gk_max(0.01*nvtxs, 25), 150); (line 240)
    // gk_max(a,b) is (a >= b ? a : b). gk_min(a,b) is (a >= b ? b : a).
    // 0.01*nvtxs is f64. 25 is int, promoted to f64. Result f64.
    // gk_max(0.01*nvtxs, 25.0) -> f64.
    // gk_min(that, 150.0) -> f64.
    // Assigned to idx_t (i32), truncating.
    let tmp_max: f64 = if (0.01 * (nvtxs as f64)) >= 25.0 {
        0.01 * (nvtxs as f64)
    } else {
        25.0
    };
    let limit: idx_t = (if tmp_max >= 150.0 { 150.0 } else { tmp_max }) as idx_t;

    // ffactor = .5/gk_max(20, nvtxs); (line 245)
    // gk_max(20, nvtxs) is int ternary. Result int.
    // .5 is double. int promoted to double. Division in double.
    // Stored in real_t (f32).
    let gk_max_20_nvtxs: idx_t = if 20 >= nvtxs { 20 } else { nvtxs };
    let ffactor: real_t = (0.5 / (gk_max_20_nvtxs as f64)) as real_t;

    // Initialize the queues (lines 248-252)
    // queues = (rpq_t **)wspacemalloc(ctrl, 2*ncon*sizeof(rpq_t *));
    // In Rust, we allocate a Vec of rpq_t directly.
    let mut queues: Vec<rpq_t> = Vec::with_capacity((2 * ncon) as usize);
    for _i in 0..(2 * ncon) {
        queues.push(super::pqueue::rpqCreate(nvtxs as usize));
    }

    for i in 0..nvtxs {
        // qnum[i] = iargmax_nrm(ncon, vwgt+i*ncon, invtvwgt); (line 252)
        let offset = (i * ncon) as usize;
        qnum[i as usize] = super::mcutil::iargmax_nrm(
            ncon as usize,
            &graph.vwgt[offset..offset + ncon as usize],
            &graph.invtvwgt,
        );
    }

    // Determine the unbalance tolerance for each constraint (lines 254-262)
    // origbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl->pijbm, ctrl->ubfactors, ubfactors); (line 260)
    // PORT NOTE: origbal is only used by the dropped debug printer (lines 265-266).
    // The call writes to ubfactors, so it must be kept.
    let _origbal: real_t = super::mcutil::ComputeLoadImbalanceDiffVec(
        graph,
        2,
        &ctrl.pijbm,
        &ctrl.ubfactors,
        &mut ubfactors,
    );

    // for (i=0; i<ncon; i++) (lines 261-262)
    for i in 0..ncon {
        // ubfactors[i] = (ubfactors[i] > 0 ? ctrl->ubfactors[i]+ubfactors[i] : ctrl->ubfactors[i]);
        if ubfactors[i as usize] > 0.0 {
            ubfactors[i as usize] = ctrl.ubfactors[i as usize] + ubfactors[i as usize];
        } else {
            ubfactors[i as usize] = ctrl.ubfactors[i as usize];
        }
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, Print2WayRefineStats(...)); (lines 265-266)
    // Dropped: debug-only under NDEBUG.

    // iset(nvtxs, -1, moved); (line 268)
    super::gklib::iset(nvtxs as usize, -1, &mut moved);

    // for (pass=0; pass<niter; pass++) (line 269)
    let mut pass: idx_t = 0;
    while pass < niter {
        // for (i=0; i<2*ncon; i++) rpqReset(queues[i]); (lines 270-271)
        for i in 0..(2 * ncon) {
            super::pqueue::rpqReset(&mut queues[i as usize]);
        }

        let mut mincutorder: idx_t = -1;
        let mut newcut: idx_t = graph.mincut;
        let mut mincut: idx_t = graph.mincut;
        let initcut: idx_t = graph.mincut;

        // minbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl->pijbm, ubfactors, minbalv); (line 276)
        let mut minbal: real_t = super::mcutil::ComputeLoadImbalanceDiffVec(
            graph,
            2,
            &ctrl.pijbm,
            &ubfactors,
            &mut minbalv,
        );

        // Insert boundary nodes in the priority queues (lines 279-288)
        let mut nbnd: idx_t = graph.nbnd;
        super::gklib::irandArrayPermute(nbnd, &mut perm, nbnd / 5, 1, &mut ctrl.rng_state);

        for ii in 0..nbnd {
            let i: idx_t = graph.bndind[perm[ii as usize] as usize];
            // rgain = ed[i]-id[i]; (line 286)
            let rgain: real_t = (graph.ed[i as usize] - graph.id[i as usize]) as real_t;
            // rpqInsert(queues[2*qnum[i]+where[i]], i, rgain); (line 287)
            let qidx: usize = (2 * qnum[i as usize] + graph.where_[i as usize]) as usize;
            super::pqueue::rpqInsert(&mut queues[qidx], i, rgain);
        }

        // for (nswaps=0; nswaps<nvtxs; nswaps++) (line 290)
        let mut nswaps: idx_t = 0;
        while nswaps < nvtxs {
            let mut from: idx_t = -1;
            let mut cnum: idx_t = -1;

            // SelectQueue(graph, ctrl->pijbm, ubfactors, queues, &from, &cnum); (line 291)
            SelectQueue(
                graph,
                &ctrl.pijbm,
                &ubfactors,
                &queues,
                &mut from,
                &mut cnum,
            );

            let to: idx_t = (from + 1) % 2;

            // if (from == -1 || (higain = rpqGetTop(queues[2*cnum+from])) == -1) break; (lines 295-296)
            let higain: idx_t;
            if from == -1 {
                break;
            } else {
                let qidx: usize = (2 * cnum + from) as usize;
                higain = super::pqueue::rpqGetTop(&mut queues[qidx]);
                if higain == -1 {
                    break;
                }
            }

            // newcut -= (ed[higain]-id[higain]); (line 298)
            newcut -= graph.ed[higain as usize] - graph.id[higain as usize];

            // iaxpy(ncon,  1, vwgt+higain*ncon, 1, pwgts+to*ncon,   1); (line 300)
            let vwgt_off = (higain * ncon) as usize;
            let pwgts_to_off = (to * ncon) as usize;
            let pwgts_from_off = (from * ncon) as usize;

            super::gklib::iaxpy(
                ncon as usize,
                1,
                &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                1,
                &mut graph.pwgts[pwgts_to_off..pwgts_to_off + ncon as usize],
                1,
            );
            // iaxpy(ncon, -1, vwgt+higain*ncon, 1, pwgts+from*ncon, 1); (line 301)
            super::gklib::iaxpy(
                ncon as usize,
                -1,
                &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                1,
                &mut graph.pwgts[pwgts_from_off..pwgts_from_off + ncon as usize],
                1,
            );

            // newbal = ComputeLoadImbalanceDiffVec(graph, 2, ctrl->pijbm, ubfactors, newbalv); (line 302)
            let newbal: real_t = super::mcutil::ComputeLoadImbalanceDiffVec(
                graph,
                2,
                &ctrl.pijbm,
                &ubfactors,
                &mut newbalv,
            );

            // if ((newcut < mincut && newbal <= ffactor) ||
            //     (newcut == mincut && (newbal < minbal ||
            //      (newbal == minbal && BetterBalance2Way(ncon, minbalv, newbalv))))) { (lines 304-306)
            let cond1 = newcut < mincut && newbal <= ffactor;
            let cond2_inner = newbal < minbal
                || (newbal == minbal
                    && super::mcutil::BetterBalance2Way(ncon, &minbalv, &newbalv) != 0);
            let cond2 = newcut == mincut && cond2_inner;

            if cond1 || cond2 {
                mincut = newcut;
                minbal = newbal;
                mincutorder = nswaps;
                // rcopy(ncon, newbalv, minbalv); (line 310)
                super::gklib::rcopy(ncon as usize, &newbalv, &mut minbalv);
            } else if nswaps - mincutorder > limit {
                // We hit the limit, undo last move (lines 312-317)
                newcut += graph.ed[higain as usize] - graph.id[higain as usize];
                // PORT NOTE: C leaves this store dead (the loop breaks and newcut is reset each pass);
                // read it once so the Rust unused_assignments lint stays quiet. No semantic effect.
                let _ = newcut;
                // iaxpy(ncon,  1, vwgt+higain*ncon, 1, pwgts+from*ncon, 1); (line 314)
                super::gklib::iaxpy(
                    ncon as usize,
                    1,
                    &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                    1,
                    &mut graph.pwgts[pwgts_from_off..pwgts_from_off + ncon as usize],
                    1,
                );
                // iaxpy(ncon, -1, vwgt+higain*ncon, 1, pwgts+to*ncon,   1); (line 315)
                super::gklib::iaxpy(
                    ncon as usize,
                    -1,
                    &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                    1,
                    &mut graph.pwgts[pwgts_to_off..pwgts_to_off + ncon as usize],
                    1,
                );
                break;
            }

            // where[higain] = to; (line 319)
            graph.where_[higain as usize] = to;
            // moved[higain] = nswaps; (line 320)
            moved[higain as usize] = nswaps;
            // swaps[nswaps] = higain; (line 321)
            swaps[nswaps as usize] = higain;

            // if (ctrl->dbglvl&METIS_DBG_MOVEINFO) { ... } (lines 323-330)
            // Dropped: debug-only under NDEBUG.

            // Update the id[i]/ed[i] values of the affected nodes (lines 333-375)
            // SWAP(id[higain], ed[higain], tmp); (line 336)
            std::mem::swap(
                &mut graph.id[higain as usize],
                &mut graph.ed[higain as usize],
            );

            // if (ed[higain] == 0 && xadj[higain] < xadj[higain+1])
            //   BNDDelete(nbnd, bndind,  bndptr, higain); (lines 337-338)
            if graph.ed[higain as usize] == 0 && xadj[higain as usize] < xadj[(higain + 1) as usize]
            {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            }

            // for (j=xadj[higain]; j<xadj[higain+1]; j++) (line 340)
            let j_start = xadj[higain as usize];
            let j_end = xadj[(higain + 1) as usize];
            let mut j = j_start;
            while j < j_end {
                let k: idx_t = adjncy[j as usize];

                // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 343)
                let kwgt: idx_t = if to == graph.where_[k as usize] {
                    graph.adjwgt[j as usize]
                } else {
                    -graph.adjwgt[j as usize]
                };

                // INC_DEC(id[k], ed[k], kwgt); (line 344)
                // id[k] += kwgt; ed[k] -= kwgt;
                graph.id[k as usize] += kwgt;
                graph.ed[k as usize] -= kwgt;

                // Update its boundary information and queue position (lines 346-374)
                if graph.bndptr[k as usize] != -1 {
                    // If k was a boundary vertex
                    if graph.ed[k as usize] == 0 {
                        // Not a boundary vertex any more
                        super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                        if moved[k as usize] == -1 {
                            // Remove it if in the queues
                            let qidx: usize =
                                (2 * qnum[k as usize] + graph.where_[k as usize]) as usize;
                            super::pqueue::rpqDelete(&mut queues[qidx], k);
                        }
                    } else {
                        // If it has not been moved, update its position in the queue
                        if moved[k as usize] == -1 {
                            // rgain = ed[k]-id[k]; (line 358)
                            let rgain: real_t =
                                (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                            let qidx: usize =
                                (2 * qnum[k as usize] + graph.where_[k as usize]) as usize;
                            super::pqueue::rpqUpdate(&mut queues[qidx], k, rgain);
                        }
                    }
                } else {
                    if graph.ed[k as usize] > 0 {
                        // It will now become a boundary vertex
                        super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                        if moved[k as usize] == -1 {
                            // rgain = ed[k]-id[k]; (line 370)
                            let rgain: real_t =
                                (graph.ed[k as usize] - graph.id[k as usize]) as real_t;
                            let qidx: usize =
                                (2 * qnum[k as usize] + graph.where_[k as usize]) as usize;
                            super::pqueue::rpqInsert(&mut queues[qidx], k, rgain);
                        }
                    }
                }
                j += 1;
            }

            nswaps += 1;
        }

        // Roll back computations (lines 380-408)
        // for (i=0; i<nswaps; i++) moved[swaps[i]] = -1; (lines 383-384)
        for i in 0..nswaps {
            moved[swaps[i as usize] as usize] = -1;
        }

        // for (nswaps--; nswaps>mincutorder; nswaps--) (line 385)
        nswaps -= 1;
        while nswaps > mincutorder {
            let higain: idx_t = swaps[nswaps as usize];

            // to = where[higain] = (where[higain]+1)%2; (line 388)
            let to: idx_t = (graph.where_[higain as usize] + 1) % 2;
            graph.where_[higain as usize] = to;

            // SWAP(id[higain], ed[higain], tmp); (line 389)
            std::mem::swap(
                &mut graph.id[higain as usize],
                &mut graph.ed[higain as usize],
            );

            // if (ed[higain] == 0 && bndptr[higain] != -1 && xadj[higain] < xadj[higain+1])
            //   BNDDelete(nbnd, bndind,  bndptr, higain); (lines 390-391)
            if graph.ed[higain as usize] == 0
                && graph.bndptr[higain as usize] != -1
                && xadj[higain as usize] < xadj[(higain + 1) as usize]
            {
                super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            } else if graph.ed[higain as usize] > 0 && graph.bndptr[higain as usize] == -1 {
                // BNDInsert(nbnd, bndind,  bndptr, higain); (lines 392-393)
                super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            }

            // iaxpy(ncon,  1, vwgt+higain*ncon, 1, pwgts+to*ncon,         1); (line 395)
            let vwgt_off = (higain * ncon) as usize;
            let pwgts_to_off = (to * ncon) as usize;
            let pwgts_other_off = (((to + 1) % 2) * ncon) as usize;

            super::gklib::iaxpy(
                ncon as usize,
                1,
                &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                1,
                &mut graph.pwgts[pwgts_to_off..pwgts_to_off + ncon as usize],
                1,
            );
            // iaxpy(ncon, -1, vwgt+higain*ncon, 1, pwgts+((to+1)%2)*ncon, 1); (line 396)
            super::gklib::iaxpy(
                ncon as usize,
                -1,
                &graph.vwgt[vwgt_off..vwgt_off + ncon as usize],
                1,
                &mut graph.pwgts[pwgts_other_off..pwgts_other_off + ncon as usize],
                1,
            );

            // for (j=xadj[higain]; j<xadj[higain+1]; j++) (line 397)
            let j_start = xadj[higain as usize];
            let j_end = xadj[(higain + 1) as usize];
            let mut j = j_start;
            while j < j_end {
                let k: idx_t = adjncy[j as usize];

                // kwgt = (to == where[k] ? adjwgt[j] : -adjwgt[j]); (line 400)
                let kwgt: idx_t = if to == graph.where_[k as usize] {
                    graph.adjwgt[j as usize]
                } else {
                    -graph.adjwgt[j as usize]
                };

                // INC_DEC(id[k], ed[k], kwgt); (line 401)
                graph.id[k as usize] += kwgt;
                graph.ed[k as usize] -= kwgt;

                // if (bndptr[k] != -1 && ed[k] == 0)
                //   BNDDelete(nbnd, bndind, bndptr, k); (lines 403-404)
                if graph.bndptr[k as usize] != -1 && graph.ed[k as usize] == 0 {
                    super::BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                }
                // if (bndptr[k] == -1 && ed[k] > 0)
                //   BNDInsert(nbnd, bndind, bndptr, k); (lines 405-406)
                if graph.bndptr[k as usize] == -1 && graph.ed[k as usize] > 0 {
                    super::BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);
                }
                j += 1;
            }

            nswaps -= 1;
        }

        // graph->mincut = mincut; (line 410)
        graph.mincut = mincut;
        // graph->nbnd   = nbnd; (line 411)
        graph.nbnd = nbnd;

        // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, Print2WayRefineStats(...)); (lines 413-414)
        // Dropped: debug-only under NDEBUG.

        // if (mincutorder <= 0 || mincut == initcut) break; (lines 416-417)
        if mincutorder <= 0 || mincut == initcut {
            break;
        }

        pass += 1;
    }

    // for (i=0; i<2*ncon; i++) rpqDestroy(queues[i]); (lines 420-421)
    for queue in queues {
        super::pqueue::rpqDestroy(queue);
    }

    // WCOREPOP (line 423)
    super::wspace::wspacepop(ctrl);
}

/// HiGHS: fm.c::SelectQueue (lines 431-493)
/// `ncon > 1` only (called by `FM_Mc2WayCutRefine`); `queues` holds the
/// `2*ncon` queues the C passes as `rpq_t **`.
pub fn SelectQueue(
    graph: &graph_t,
    pijbm: &[real_t],
    ubfactors: &[real_t],
    queues: &[rpq_t],
    from: &mut idx_t,
    cnum: &mut idx_t,
) {
    let ncon = graph.ncon;

    *from = -1;
    *cnum = -1;

    // First determine the side and the queue, irrespective of the presence of nodes.
    // The side & queue is determined based on the most violated balancing constraint.
    // for (max=0.0, part=0; part<2; part++) { (line 444)
    let mut max: real_t = 0.0;
    for part in 0..2 {
        // for (i=0; i<ncon; i++) { (line 445)
        for i in 0..ncon {
            // tmp = graph->pwgts[part*ncon+i]*pijbm[part*ncon+i] - ubfactors[i]; (line 446)
            // pwgts is idx_t (i32), pijbm is real_t (f32). C: i32 * f32 -> f32.
            // Then f32 - f32 -> f32.
            let tmp: real_t = (graph.pwgts[(part * ncon + i) as usize] as f32)
                * pijbm[(part * ncon + i) as usize]
                - ubfactors[i as usize];

            // if (tmp >= max) { (line 449)
            if tmp >= max {
                max = tmp;
                *from = part;
                *cnum = i;
            }
        }
    }

    // if (*from != -1) { (line 458)
    if *from != -1 {
        // in case the desired queue is empty, select a queue from the same side
        // if (rpqLength(queues[2*(*cnum)+(*from)]) == 0) { (line 460)
        if super::pqueue::rpqLength(&queues[(2 * (*cnum) + (*from)) as usize]) == 0 {
            // for (i=0; i<ncon; i++) { (line 461)
            let mut i: idx_t = 0;
            while i < ncon {
                // if (rpqLength(queues[2*i+(*from)]) > 0) { (line 462)
                if super::pqueue::rpqLength(&queues[(2 * i + (*from)) as usize]) > 0 {
                    // max = graph->pwgts[(*from)*ncon+i]*pijbm[(*from)*ncon+i] - ubfactors[i]; (line 463)
                    max = (graph.pwgts[((*from) * ncon + i) as usize] as f32)
                        * pijbm[((*from) * ncon + i) as usize]
                        - ubfactors[i as usize];
                    *cnum = i;
                    break;
                }
                i += 1;
            }

            // for (i++; i<ncon; i++) { (line 469)
            i += 1;
            while i < ncon {
                // tmp = graph->pwgts[(*from)*ncon+i]*pijbm[(*from)*ncon+i] - ubfactors[i]; (line 470)
                let tmp: real_t = (graph.pwgts[((*from) * ncon + i) as usize] as f32)
                    * pijbm[((*from) * ncon + i) as usize]
                    - ubfactors[i as usize];

                // if (tmp > max && rpqLength(queues[2*i+(*from)]) > 0) { (line 471)
                if tmp > max && super::pqueue::rpqLength(&queues[(2 * i + (*from)) as usize]) > 0 {
                    max = tmp;
                    *cnum = i;
                }
                i += 1;
            }
        }
    } else {
        // the partitioning does not violate balancing constraints, in which case select
        // a queue based on cut criteria
        // for (part=0; part<2; part++) { (line 482)
        for part in 0..2 {
            // for (i=0; i<ncon; i++) { (line 483)
            for i in 0..ncon {
                // if (rpqLength(queues[2*i+part]) > 0 &&
                //     (*from == -1 || rpqSeeTopKey(queues[2*i+part]) > max)) { (lines 484-485)
                if super::pqueue::rpqLength(&queues[(2 * i + part) as usize]) > 0
                    && (*from == -1
                        || super::pqueue::rpqSeeTopKey(&queues[(2 * i + part) as usize]) > max)
                {
                    // max = rpqSeeTopKey(queues[2*i+part]); (line 486)
                    max = super::pqueue::rpqSeeTopKey(&queues[(2 * i + part) as usize]);
                    *from = part;
                    *cnum = i;
                }
            }
        }
    }
}

/// HiGHS: fm.c::Print2WayRefineStats (lines 499-526)
/// Debug printer (`METIS_DBG_MOVEINFO`); never reached from HiPO.
pub fn Print2WayRefineStats(
    ctrl: &ctrl_t,
    graph: &graph_t,
    ntpwgts: &[real_t],
    deltabal: real_t,
    mincutorder: idx_t,
) {
    // HIGHS_ORDERING_PRINT expands to nothing under NDEBUG (oracle build).
    // C lines 499-526: all output is debug-only.
    let _ = ctrl;
    let _ = graph;
    let _ = ntpwgts;
    let _ = deltabal;
    let _ = mincutorder;
}
