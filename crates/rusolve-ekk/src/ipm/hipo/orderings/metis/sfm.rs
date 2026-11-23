//! HiGHS source: extern/metis/libmetis/sfm.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Node-separator FM refinement (2-sided and 1-sided) and node balancing. The
//! two `rpq_t` queues and the `rng_state`-driven `irandArrayPermute` decide the
//! move order.

use super::gklib::irandArrayPermute;
use super::pqueue::{
    rpq_t, rpqCreate, rpqDelete, rpqDestroy, rpqGetTop, rpqInsert, rpqReset, rpqSeeTopVal,
    rpqUpdate,
};
use super::{BNDDelete, BNDInsert};
use super::{ctrl_t, graph_t, idx_t, real_t};

/// HiGHS: sfm.c::FM_2WayNodeRefine2Sided (lines 21-245)
pub fn FM_2WayNodeRefine2Sided(ctrl: &mut ctrl_t, graph: &mut graph_t, niter: idx_t) {
    // WCOREPUSH; // level-0 no-op, dropped

    let nvtxs = graph.nvtxs;
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();
    let vwgt: Vec<idx_t> = graph.vwgt.clone();

    let mut queues: [rpq_t; 2] = [rpqCreate(nvtxs as usize), rpqCreate(nvtxs as usize)];

    // iwspacemalloc -> fresh zeroed Vecs
    let mut moved = vec![-1 as idx_t; nvtxs as usize];
    let mut swaps = vec![0 as idx_t; nvtxs as usize];
    let mut mptr = vec![0 as idx_t; (nvtxs + 1) as usize];
    let mut mind = vec![0 as idx_t; (2 * nvtxs) as usize];

    // mult = 0.5*ctrl->ubfactors[0];
    // ubfactors[0] is real_t (f32). 0.5 is a double literal.
    // C: 0.5 * (f32) -> promotes f32 to f64, multiplies in f64, result f64.
    // mult is real_t (f32), so the f64 result is rounded to f32.
    let mult: real_t = (0.5 * (ctrl.ubfactors[0] as f64)) as real_t;

    // badmaxpwgt = (idx_t)(mult*(pwgts[0]+pwgts[1]+pwgts[2]));
    // pwgts[0]+pwgts[1]+pwgts[2] is idx_t (i32) arithmetic.
    // mult is real_t (f32). C: f32 * i32 -> f32. Result cast to idx_t.
    let sum_pwgts = graph.pwgts[0] + graph.pwgts[1] + graph.pwgts[2];
    let badmaxpwgt: idx_t = (mult * (sum_pwgts as real_t)) as idx_t;

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 59-60) dropped

    let mut pass: idx_t = 0;
    while pass < niter {
        // iset(nvtxs, -1, moved);
        for i in 0..nvtxs as usize {
            moved[i] = -1;
        }
        rpqReset(&mut queues[0]);
        rpqReset(&mut queues[1]);

        let mut mincutorder: idx_t = -1;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        let mut nbnd = graph.nbnd;

        // irandArrayPermute(nbnd, swaps, nbnd, 1, &ctrl->rng_state);
        irandArrayPermute(nbnd, &mut swaps, nbnd, 1, &mut ctrl.rng_state);

        for ii in 0..nbnd {
            let i = graph.bndind[swaps[ii as usize] as usize];
            // rpqInsert(queues[0], i, vwgt[i]-rinfo[i].edegrees[1]);
            let key0 = (vwgt[i as usize] - graph.nrinfo[i as usize].edegrees[1]) as real_t;
            rpqInsert(&mut queues[0], i, key0);
            // rpqInsert(queues[1], i, vwgt[i]-rinfo[i].edegrees[0]);
            let key1 = (vwgt[i as usize] - graph.nrinfo[i as usize].edegrees[0]) as real_t;
            rpqInsert(&mut queues[1], i, key1);
        }

        // limit = (ctrl->compress ? gk_min(5*nbnd, 400) : gk_min(2*nbnd, 300));
        // gk_min is a C ternary: ((a) >= (b) ? (b) : (a))
        let limit: idx_t = if ctrl.compress != 0 {
            let a = 5 * nbnd;
            let b = 400;
            if a >= b { b } else { a }
        } else {
            let a = 2 * nbnd;
            let b = 300;
            if a >= b { b } else { a }
        };

        mptr[0] = 0;
        let mut nmind: idx_t = 0;
        let mut mindiff: idx_t = (graph.pwgts[0] - graph.pwgts[1]).abs();
        // C: `to = (pwgts[0] < pwgts[1] ? 0 : 1)` is dead here: every path of the
        // move loop below assigns `to` before reading it.
        let mut to: idx_t;

        let mut nswaps: idx_t = 0;
        while nswaps < nvtxs {
            let u0 = rpqSeeTopVal(&queues[0]);
            let u1 = rpqSeeTopVal(&queues[1]);

            if u0 != -1 && u1 != -1 {
                let g0 = vwgt[u0 as usize] - graph.nrinfo[u0 as usize].edegrees[1];
                let g1 = vwgt[u1 as usize] - graph.nrinfo[u1 as usize].edegrees[0];

                to = if g0 > g1 {
                    0
                } else if g0 < g1 {
                    1
                } else {
                    pass % 2
                };

                // if (pwgts[to]+vwgt[u[to]] > badmaxpwgt) to = (to+1)%2;
                let u_to = if to == 0 { u0 } else { u1 };
                if graph.pwgts[to as usize] + vwgt[u_to as usize] > badmaxpwgt {
                    to = (to + 1) % 2;
                }
            } else if u0 == -1 && u1 == -1 {
                break;
            } else if u0 != -1 && graph.pwgts[0] + vwgt[u0 as usize] <= badmaxpwgt {
                to = 0;
            } else if u1 != -1 && graph.pwgts[1] + vwgt[u1 as usize] <= badmaxpwgt {
                to = 1;
            } else {
                break;
            }

            let other = (to + 1) % 2;

            let higain = rpqGetTop(&mut queues[to as usize]);
            if moved[higain as usize] == -1 {
                rpqDelete(&mut queues[other as usize], higain);
            }

            // if (nmind + xadj[higain+1]-xadj[higain] >= 2*nvtxs-1) break;
            if nmind + (xadj[(higain + 1) as usize] - xadj[higain as usize]) >= 2 * nvtxs - 1 {
                break;
            }

            // pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);
            graph.pwgts[2] -=
                vwgt[higain as usize] - graph.nrinfo[higain as usize].edegrees[other as usize];

            // newdiff = iabs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
            let newdiff: idx_t = (graph.pwgts[to as usize] + vwgt[higain as usize]
                - (graph.pwgts[other as usize]
                    - graph.nrinfo[higain as usize].edegrees[other as usize]))
                .abs();

            if graph.pwgts[2] < mincut || (graph.pwgts[2] == mincut && newdiff < mindiff) {
                mincut = graph.pwgts[2];
                mincutorder = nswaps;
                mindiff = newdiff;
            } else {
                // if (nswaps - mincutorder > 2*limit || (nswaps - mincutorder > limit && pwgts[2] > 1.10*mincut))
                if nswaps - mincutorder > 2 * limit
                    || (nswaps - mincutorder > limit
                        && (graph.pwgts[2] as f64) > 1.10 * (mincut as f64))
                {
                    graph.pwgts[2] += vwgt[higain as usize]
                        - graph.nrinfo[higain as usize].edegrees[other as usize];
                    break;
                }
            }

            BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            graph.pwgts[to as usize] += vwgt[higain as usize];
            graph.where_[higain as usize] = to;
            moved[higain as usize] = nswaps;
            swaps[nswaps as usize] = higain;

            // Update the degrees of the affected nodes
            let mut j = xadj[higain as usize];
            while j < xadj[(higain + 1) as usize] {
                let k = adjncy[j as usize];
                if graph.where_[k as usize] == 2 {
                    // For the in-separator vertices modify their edegree[to]
                    let oldgain = vwgt[k as usize] - graph.nrinfo[k as usize].edegrees[to as usize];
                    graph.nrinfo[k as usize].edegrees[to as usize] += vwgt[higain as usize];
                    if moved[k as usize] == -1 || moved[k as usize] == -(2 + other) {
                        let newkey = (oldgain - vwgt[higain as usize]) as real_t;
                        rpqUpdate(&mut queues[other as usize], k, newkey);
                    }
                } else if graph.where_[k as usize] == other {
                    // This vertex is pulled into the separator
                    BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);

                    mind[nmind as usize] = k;
                    nmind += 1;
                    graph.where_[k as usize] = 2;
                    graph.pwgts[other as usize] -= vwgt[k as usize];

                    graph.nrinfo[k as usize].edegrees[0] = 0;
                    graph.nrinfo[k as usize].edegrees[1] = 0;

                    let mut jj = xadj[k as usize];
                    while jj < xadj[(k + 1) as usize] {
                        let kk = adjncy[jj as usize];
                        if graph.where_[kk as usize] != 2 {
                            graph.nrinfo[k as usize].edegrees
                                [graph.where_[kk as usize] as usize] += vwgt[kk as usize];
                        } else {
                            let oldgain = vwgt[kk as usize]
                                - graph.nrinfo[kk as usize].edegrees[other as usize];
                            graph.nrinfo[kk as usize].edegrees[other as usize] -= vwgt[k as usize];
                            if moved[kk as usize] == -1 || moved[kk as usize] == -(2 + to) {
                                let newkey = (oldgain + vwgt[k as usize]) as real_t;
                                rpqUpdate(&mut queues[to as usize], kk, newkey);
                            }
                        }
                        jj += 1;
                    }

                    // Insert the new vertex into the priority queue. Only one side!
                    if moved[k as usize] == -1 {
                        let newkey = (vwgt[k as usize]
                            - graph.nrinfo[k as usize].edegrees[other as usize])
                            as real_t;
                        rpqInsert(&mut queues[to as usize], k, newkey);
                        moved[k as usize] = -(2 + to);
                    }
                }
                j += 1;
            }
            mptr[(nswaps + 1) as usize] = nmind;

            // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...) (sfm.c 188-189) dropped

            nswaps += 1;
        }

        // Roll back computation
        nswaps -= 1;
        while nswaps > mincutorder {
            let higain = swaps[nswaps as usize];

            to = graph.where_[higain as usize];
            let other = (to + 1) % 2;

            // INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
            graph.pwgts[2] += vwgt[higain as usize];
            graph.pwgts[to as usize] -= vwgt[higain as usize];

            graph.where_[higain as usize] = 2;
            BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);

            graph.nrinfo[higain as usize].edegrees[0] = 0;
            graph.nrinfo[higain as usize].edegrees[1] = 0;

            let mut j = xadj[higain as usize];
            while j < xadj[(higain + 1) as usize] {
                let k = adjncy[j as usize];
                if graph.where_[k as usize] == 2 {
                    graph.nrinfo[k as usize].edegrees[to as usize] -= vwgt[higain as usize];
                } else {
                    graph.nrinfo[higain as usize].edegrees[graph.where_[k as usize] as usize] +=
                        vwgt[k as usize];
                }
                j += 1;
            }

            // Push nodes out of the separator
            let mut j = mptr[nswaps as usize];
            while j < mptr[(nswaps + 1) as usize] {
                let k = mind[j as usize];
                graph.where_[k as usize] = other;
                // INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
                graph.pwgts[other as usize] += vwgt[k as usize];
                graph.pwgts[2] -= vwgt[k as usize];
                BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);

                let mut jj = xadj[k as usize];
                while jj < xadj[(k + 1) as usize] {
                    let kk = adjncy[jj as usize];
                    if graph.where_[kk as usize] == 2 {
                        graph.nrinfo[kk as usize].edegrees[other as usize] += vwgt[k as usize];
                    }
                    jj += 1;
                }
                j += 1;
            }
            nswaps -= 1;
        }

        // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 231-232) dropped

        graph.mincut = mincut;
        graph.nbnd = nbnd;

        if mincutorder == -1 || mincut >= initcut {
            break;
        }

        pass += 1;
    }

    let [q0, q1] = queues;
    rpqDestroy(q0);
    rpqDestroy(q1);

    // WCOREPOP; // level-0 no-op, dropped
}

/// HiGHS: sfm.c::FM_2WayNodeRefine1Sided (lines 255-444)
pub fn FM_2WayNodeRefine1Sided(ctrl: &mut ctrl_t, graph: &mut graph_t, niter: idx_t) {
    // WCOREPUSH; // level-0 no-op, dropped

    let nvtxs = graph.nvtxs;
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();
    let vwgt: Vec<idx_t> = graph.vwgt.clone();

    let mut queue = rpqCreate(nvtxs as usize);

    // iwspacemalloc -> fresh zeroed Vecs
    let mut swaps = vec![0 as idx_t; nvtxs as usize];
    let mut mptr = vec![0 as idx_t; (nvtxs + 1) as usize];
    let mut mind = vec![0 as idx_t; (2 * nvtxs) as usize];

    // mult = 0.5*ctrl->ubfactors[0];
    let mult: real_t = (0.5 * (ctrl.ubfactors[0] as f64)) as real_t;

    // badmaxpwgt = (idx_t)(mult*(pwgts[0]+pwgts[1]+pwgts[2]));
    let sum_pwgts = graph.pwgts[0] + graph.pwgts[1] + graph.pwgts[2];
    let badmaxpwgt: idx_t = (mult * (sum_pwgts as real_t)) as idx_t;

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 290-291) dropped

    let mut to: idx_t = if graph.pwgts[0] < graph.pwgts[1] {
        1
    } else {
        0
    };

    let mut pass: idx_t = 0;
    while pass < 2 * niter {
        let other = to;
        to = (to + 1) % 2;

        rpqReset(&mut queue);

        let mut mincutorder: idx_t = -1;
        let initcut = graph.mincut;
        let mut mincut = graph.mincut;
        let mut nbnd = graph.nbnd;

        // irandArrayPermute(nbnd, swaps, nbnd, 1, &ctrl->rng_state);
        irandArrayPermute(nbnd, &mut swaps, nbnd, 1, &mut ctrl.rng_state);

        for ii in 0..nbnd {
            let i = graph.bndind[swaps[ii as usize] as usize];
            // rpqInsert(queue, i, vwgt[i]-rinfo[i].edegrees[other]);
            let key =
                (vwgt[i as usize] - graph.nrinfo[i as usize].edegrees[other as usize]) as real_t;
            rpqInsert(&mut queue, i, key);
        }

        // limit = (ctrl->compress ? gk_min(5*nbnd, 500) : gk_min(3*nbnd, 300));
        let limit: idx_t = if ctrl.compress != 0 {
            let a = 5 * nbnd;
            let b = 500;
            if a >= b { b } else { a }
        } else {
            let a = 3 * nbnd;
            let b = 300;
            if a >= b { b } else { a }
        };

        mptr[0] = 0;
        let mut nmind: idx_t = 0;
        let mut mindiff: idx_t = (graph.pwgts[0] - graph.pwgts[1]).abs();

        let mut nswaps: idx_t = 0;
        while nswaps < nvtxs {
            let higain = rpqGetTop(&mut queue);
            if higain == -1 {
                break;
            }

            // if (nmind + xadj[higain+1]-xadj[higain] >= 2*nvtxs-1) break;
            if nmind + (xadj[(higain + 1) as usize] - xadj[higain as usize]) >= 2 * nvtxs - 1 {
                break;
            }

            // if (pwgts[to]+vwgt[higain] > badmaxpwgt) break;
            if graph.pwgts[to as usize] + vwgt[higain as usize] > badmaxpwgt {
                break;
            }

            // pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);
            graph.pwgts[2] -=
                vwgt[higain as usize] - graph.nrinfo[higain as usize].edegrees[other as usize];

            // newdiff = iabs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
            let newdiff: idx_t = (graph.pwgts[to as usize] + vwgt[higain as usize]
                - (graph.pwgts[other as usize]
                    - graph.nrinfo[higain as usize].edegrees[other as usize]))
                .abs();

            if graph.pwgts[2] < mincut || (graph.pwgts[2] == mincut && newdiff < mindiff) {
                mincut = graph.pwgts[2];
                mincutorder = nswaps;
                mindiff = newdiff;
            } else {
                // if (nswaps - mincutorder > 3*limit || (nswaps - mincutorder > limit && pwgts[2] > 1.10*mincut))
                if nswaps - mincutorder > 3 * limit
                    || (nswaps - mincutorder > limit
                        && (graph.pwgts[2] as f64) > 1.10 * (mincut as f64))
                {
                    graph.pwgts[2] += vwgt[higain as usize]
                        - graph.nrinfo[higain as usize].edegrees[other as usize];
                    break;
                }
            }

            BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
            graph.pwgts[to as usize] += vwgt[higain as usize];
            graph.where_[higain as usize] = to;
            swaps[nswaps as usize] = higain;

            // Update the degrees of the affected nodes
            let mut j = xadj[higain as usize];
            while j < xadj[(higain + 1) as usize] {
                let k = adjncy[j as usize];

                if graph.where_[k as usize] == 2 {
                    // For the in-separator vertices modify their edegree[to]
                    graph.nrinfo[k as usize].edegrees[to as usize] += vwgt[higain as usize];
                } else if graph.where_[k as usize] == other {
                    // This vertex is pulled into the separator
                    BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);

                    mind[nmind as usize] = k;
                    nmind += 1;
                    graph.where_[k as usize] = 2;
                    graph.pwgts[other as usize] -= vwgt[k as usize];

                    graph.nrinfo[k as usize].edegrees[0] = 0;
                    graph.nrinfo[k as usize].edegrees[1] = 0;

                    let mut jj = xadj[k as usize];
                    let iend = xadj[(k + 1) as usize];
                    while jj < iend {
                        let kk = adjncy[jj as usize];
                        if graph.where_[kk as usize] != 2 {
                            graph.nrinfo[k as usize].edegrees
                                [graph.where_[kk as usize] as usize] += vwgt[kk as usize];
                        } else {
                            graph.nrinfo[kk as usize].edegrees[other as usize] -= vwgt[k as usize];

                            // rpqUpdate(queue, kk, vwgt[kk]-rinfo[kk].edegrees[other]);
                            let newkey = (vwgt[kk as usize]
                                - graph.nrinfo[kk as usize].edegrees[other as usize])
                                as real_t;
                            rpqUpdate(&mut queue, kk, newkey);
                        }
                        jj += 1;
                    }

                    // Insert the new vertex into the priority queue. Safe due to one-sided moves
                    let newkey = (vwgt[k as usize]
                        - graph.nrinfo[k as usize].edegrees[other as usize])
                        as real_t;
                    rpqInsert(&mut queue, k, newkey);
                }
                j += 1;
            }
            mptr[(nswaps + 1) as usize] = nmind;

            // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...) (sfm.c 390-393) dropped

            nswaps += 1;
        }

        // Roll back computation
        nswaps -= 1;
        while nswaps > mincutorder {
            let higain = swaps[nswaps as usize];

            // INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
            graph.pwgts[2] += vwgt[higain as usize];
            graph.pwgts[to as usize] -= vwgt[higain as usize];

            graph.where_[higain as usize] = 2;
            BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);

            graph.nrinfo[higain as usize].edegrees[0] = 0;
            graph.nrinfo[higain as usize].edegrees[1] = 0;

            let mut j = xadj[higain as usize];
            while j < xadj[(higain + 1) as usize] {
                let k = adjncy[j as usize];
                if graph.where_[k as usize] == 2 {
                    graph.nrinfo[k as usize].edegrees[to as usize] -= vwgt[higain as usize];
                } else {
                    graph.nrinfo[higain as usize].edegrees[graph.where_[k as usize] as usize] +=
                        vwgt[k as usize];
                }
                j += 1;
            }

            // Push nodes out of the separator
            let mut j = mptr[nswaps as usize];
            while j < mptr[(nswaps + 1) as usize] {
                let k = mind[j as usize];
                graph.where_[k as usize] = other;
                // INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
                graph.pwgts[other as usize] += vwgt[k as usize];
                graph.pwgts[2] -= vwgt[k as usize];
                BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);

                let mut jj = xadj[k as usize];
                let iend = xadj[(k + 1) as usize];
                while jj < iend {
                    let kk = adjncy[jj as usize];
                    if graph.where_[kk as usize] == 2 {
                        graph.nrinfo[kk as usize].edegrees[other as usize] += vwgt[k as usize];
                    }
                    jj += 1;
                }
                j += 1;
            }
            nswaps -= 1;
        }

        // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 431-432) dropped

        graph.mincut = mincut;
        graph.nbnd = nbnd;

        if pass % 2 == 1 && (mincutorder == -1 || mincut >= initcut) {
            break;
        }

        pass += 1;
    }

    rpqDestroy(queue);

    // WCOREPOP; // level-0 no-op, dropped
}

/// HiGHS: sfm.c::FM_2WayNodeBalance (lines 451-580)
pub fn FM_2WayNodeBalance(ctrl: &mut ctrl_t, graph: &mut graph_t) {
    // WCOREPUSH; // level-0 no-op, dropped

    let nvtxs = graph.nvtxs;
    let xadj: Vec<idx_t> = graph.xadj_eff().to_vec();
    let adjncy: Vec<idx_t> = graph.adjncy_eff().to_vec();
    let vwgt: Vec<idx_t> = graph.vwgt.clone();

    // mult = 0.5*ctrl->ubfactors[0];
    // ubfactors[0] is real_t (f32). 0.5 is a double literal.
    // C: 0.5 * (f32) -> promotes f32 to f64, multiplies in f64, result f64.
    // mult is real_t (f32), so the f64 result is rounded to f32.
    let mult: real_t = (0.5 * (ctrl.ubfactors[0] as f64)) as real_t;

    // badmaxpwgt = (idx_t)(mult*(pwgts[0]+pwgts[1]));
    // pwgts[0]+pwgts[1] is idx_t (i32) arithmetic.
    // mult is real_t (f32). C: f32 * i32 -> f32. Result cast to idx_t.
    let sum_pwgts = graph.pwgts[0] + graph.pwgts[1];
    let mut badmaxpwgt: idx_t = (mult * (sum_pwgts as real_t)) as idx_t;

    // C: if (gk_max(pwgts[0], pwgts[1]) < badmaxpwgt) return;
    // gk_max is a C ternary: ((a) > (b) ? (a) : (b))
    let gk_max_pw = if graph.pwgts[0] > graph.pwgts[1] {
        graph.pwgts[0]
    } else {
        graph.pwgts[1]
    };
    if gk_max_pw < badmaxpwgt {
        return;
    }

    // C: if (iabs(pwgts[0]-pwgts[1]) < 3*graph->tvwgt[0]/nvtxs) return;
    // 3*graph->tvwgt[0] is idx_t * idx_t -> idx_t. / nvtxs is idx_t / idx_t -> idx_t (trunc).
    // iabs(...) is idx_t. Comparison idx_t < idx_t.
    if (graph.pwgts[0] - graph.pwgts[1]).abs() < 3 * graph.tvwgt[0] / nvtxs {
        return;
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 491-492) dropped

    // to = (pwgts[0] < pwgts[1] ? 0 : 1);
    let to: idx_t = if graph.pwgts[0] < graph.pwgts[1] {
        0
    } else {
        1
    };
    let other: idx_t = (to + 1) % 2;

    let mut queue = rpqCreate(nvtxs as usize);

    // iwspacemalloc -> fresh Vecs
    let mut perm = vec![0 as idx_t; nvtxs as usize];
    let mut moved = vec![-1 as idx_t; nvtxs as usize];

    let mut nbnd = graph.nbnd;
    // irandArrayPermute(nbnd, perm, nbnd, 1, &ctrl->rng_state);
    irandArrayPermute(nbnd, &mut perm, nbnd, 1, &mut ctrl.rng_state);
    for ii in 0..nbnd {
        let i = graph.bndind[perm[ii as usize] as usize];
        // rpqInsert(queue, i, vwgt[i]-rinfo[i].edegrees[other]);
        let key = (vwgt[i as usize] - graph.nrinfo[i as usize].edegrees[other as usize]) as real_t;
        rpqInsert(&mut queue, i, key);
    }

    /******************************************************
     * Get into the FM loop
     *******************************************************/
    let mut nswaps: idx_t = 0;
    while nswaps < nvtxs {
        let higain = rpqGetTop(&mut queue);
        if higain == -1 {
            break;
        }

        moved[higain as usize] = 1;

        let gain = vwgt[higain as usize] - graph.nrinfo[higain as usize].edegrees[other as usize];
        // badmaxpwgt = (idx_t)(mult*(pwgts[0]+pwgts[1]));
        let sum_pwgts = graph.pwgts[0] + graph.pwgts[1];
        badmaxpwgt = (mult * (sum_pwgts as real_t)) as idx_t;

        /* break if other is now underwight */
        if graph.pwgts[to as usize] > graph.pwgts[other as usize] {
            break;
        }

        /* break if balance is achieved and no +ve or zero gain */
        if gain < 0 && graph.pwgts[other as usize] < badmaxpwgt {
            break;
        }

        /* skip this vertex if it will violate balance on the other side */
        if graph.pwgts[to as usize] + vwgt[higain as usize] > badmaxpwgt {
            nswaps += 1;
            continue;
        }

        graph.pwgts[2] -= gain;

        BNDDelete(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, higain);
        graph.pwgts[to as usize] += vwgt[higain as usize];
        graph.where_[higain as usize] = to;

        // IFSET(ctrl->dbglvl, METIS_DBG_MOVEINFO, ...) (sfm.c 532-533) dropped

        /**********************************************************
         * Update the degrees of the affected nodes
         ***********************************************************/
        let mut j = xadj[higain as usize];
        while j < xadj[(higain + 1) as usize] {
            let k = adjncy[j as usize];
            if graph.where_[k as usize] == 2 {
                /* For the in-separator vertices modify their edegree[to] */
                graph.nrinfo[k as usize].edegrees[to as usize] += vwgt[higain as usize];
            } else if graph.where_[k as usize] == other {
                /* This vertex is pulled into the separator */
                BNDInsert(&mut nbnd, &mut graph.bndind, &mut graph.bndptr, k);

                graph.where_[k as usize] = 2;
                graph.pwgts[other as usize] -= vwgt[k as usize];

                graph.nrinfo[k as usize].edegrees[0] = 0;
                graph.nrinfo[k as usize].edegrees[1] = 0;
                let mut jj = xadj[k as usize];
                while jj < xadj[(k + 1) as usize] {
                    let kk = adjncy[jj as usize];
                    if graph.where_[kk as usize] != 2 {
                        graph.nrinfo[k as usize].edegrees[graph.where_[kk as usize] as usize] +=
                            vwgt[kk as usize];
                    } else {
                        let oldgain =
                            vwgt[kk as usize] - graph.nrinfo[kk as usize].edegrees[other as usize];
                        graph.nrinfo[kk as usize].edegrees[other as usize] -= vwgt[k as usize];

                        if moved[kk as usize] == -1 {
                            // rpqUpdate(queue, kk, oldgain+vwgt[k]);
                            let newkey = (oldgain + vwgt[k as usize]) as real_t;
                            rpqUpdate(&mut queue, kk, newkey);
                        }
                    }
                    jj += 1;
                }

                /* Insert the new vertex into the priority queue */
                // rpqInsert(queue, k, vwgt[k]-edegrees[other]);
                let newkey = (vwgt[k as usize] - graph.nrinfo[k as usize].edegrees[other as usize])
                    as real_t;
                rpqInsert(&mut queue, k, newkey);
            }
            j += 1;
        }

        nswaps += 1;
    }

    // IFSET(ctrl->dbglvl, METIS_DBG_REFINE, ...) (sfm.c 571-572) dropped

    graph.mincut = graph.pwgts[2];
    graph.nbnd = nbnd;

    rpqDestroy(queue);

    // WCOREPOP; // level-0 no-op, dropped
}
