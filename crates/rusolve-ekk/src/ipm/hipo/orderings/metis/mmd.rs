//! HiGHS source: extern/metis/libmetis/mmd.c
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! Multiple minimum degree (Liu's Fortran `genmmd` translated to C). The C
//! shifts every array pointer by one (`xadj--` ...) and indexes from 1; the
//! Rust functions take the unshifted slices and index `a[(i - 1) as usize]`
//! through one local accessor, never by allocating a padded copy.
//! `mmdelm` / `mmdupd` / `mmdint` / `mmdnum`
//! receive the already-shifted pointers in C; in Rust they take the same
//! unshifted slices.

use super::idx_t;

/// HiGHS: mmd.c::genmmd (lines 53-153)
#[allow(unused_assignments)] // C locals are initialised before use on every goto path
pub fn genmmd(
    neqns: idx_t,
    xadj: &mut [idx_t],
    adjncy: &mut [idx_t],
    invp: &mut [idx_t],
    perm: &mut [idx_t],
    delta: idx_t,
    head: &mut [idx_t],
    qsize: &mut [idx_t],
    list: &mut [idx_t],
    marker: &mut [idx_t],
    maxint: idx_t,
    ncsub: &mut idx_t,
) {
    let mut ehead: idx_t = 0;
    let mut i: idx_t = 0;
    let mut mdeg: idx_t = 0;
    let mut mdlmt: idx_t = 0;
    let mut mdeg_node: idx_t = 0;
    let mut nextmd: idx_t = 0;
    let mut num: idx_t = 0;
    let mut tag: idx_t = 0;

    if neqns <= 0 {
        return;
    }

    /* adjust from C to Fortran */
    // C: xadj--; adjncy--; invp--; perm--; head--; qsize--; list--; marker--;
    // In Rust we keep the unshifted slices and index a[(i - 1) as usize].

    /* initialization for the minimum degree algorithm */
    *ncsub = 0;
    mmdint(neqns, xadj, adjncy, head, invp, perm, qsize, list, marker);

    /* 'num' counts the number of ordered nodes plus 1 */
    num = 1;

    /* eliminate all isolated nodes */
    nextmd = head[0];
    while nextmd > 0 {
        mdeg_node = nextmd;
        nextmd = invp[(mdeg_node - 1) as usize];
        marker[(mdeg_node - 1) as usize] = maxint;
        invp[(mdeg_node - 1) as usize] = -num;
        num += 1;
    }

    /* search for node of the minimum degree. 'mdeg' is the current */
    /* minimum degree; 'tag' is used to facilitate marking nodes.   */
    if num > neqns {
        // goto n1000
        mmdnum(neqns, perm, invp, qsize);
        return;
    }
    tag = 1;
    head[0] = 0;
    mdeg = 2;

    // State machine for goto labels:
    // 0 = top of infinite loop (line 91)
    // 1 = n500 (line 102)
    // 2 = n900 (line 140)
    // 3 = n1000 (line 148)

    let mut state: i32 = 0;

    'goto: loop {
        match state {
            0 => {
                // while (1) {  line 91
                // while (head[mdeg] <= 0) mdeg++;  lines 92-93
                while head[(mdeg - 1) as usize] <= 0 {
                    mdeg += 1;
                }

                /* use value of 'delta' to set up 'mdlmt', which governs */
                /* when a degree update is to be performed.              */
                // mdlmt = mdeg + delta;  (commented out in C, line 97)
                // the need for gk_min() was identified by jsf67
                // gk_min(a, b) = ((a) >= (b) ? (b) : (a))
                mdlmt = if neqns >= mdeg + delta {
                    mdeg + delta
                } else {
                    neqns
                };
                ehead = 0;

                state = 1; // goto n500
            }
            1 => {
                // n500:  line 102
                mdeg_node = head[(mdeg - 1) as usize];
                while mdeg_node <= 0 {
                    mdeg += 1;

                    if mdeg > mdlmt {
                        state = 2; // goto n900
                        continue 'goto;
                    }
                    mdeg_node = head[(mdeg - 1) as usize];
                }

                /* remove 'mdeg_node' from the degree structure */
                nextmd = invp[(mdeg_node - 1) as usize];
                head[(mdeg - 1) as usize] = nextmd;
                if nextmd > 0 {
                    perm[(nextmd - 1) as usize] = -mdeg;
                }
                invp[(mdeg_node - 1) as usize] = -num;
                *ncsub += mdeg + qsize[(mdeg_node - 1) as usize] - 2;
                if (num + qsize[(mdeg_node - 1) as usize]) > neqns {
                    state = 3; // goto n1000
                    continue 'goto;
                }

                /*  eliminate 'mdeg_node' and perform quotient graph */
                /*  transformation. reset 'tag' value if necessary.    */
                tag += 1;
                if tag >= maxint {
                    tag = 1;
                    i = 1;
                    while i <= neqns {
                        if marker[(i - 1) as usize] < maxint {
                            marker[(i - 1) as usize] = 0;
                        }
                        i += 1;
                    }
                }

                mmdelm(
                    mdeg_node, xadj, adjncy, head, invp, perm, qsize, list, marker, maxint, tag,
                );

                num += qsize[(mdeg_node - 1) as usize];
                list[(mdeg_node - 1) as usize] = ehead;
                ehead = mdeg_node;
                if delta >= 0 {
                    state = 1; // goto n500
                    continue 'goto;
                }

                state = 2; // fall through to n900
            }
            2 => {
                // n900:  line 140
                /* update degrees of the nodes involved in the  */
                /* minimum degree nodes elimination.            */
                if num > neqns {
                    state = 3; // goto n1000
                    continue 'goto;
                }
                mmdupd(
                    ehead, neqns, xadj, adjncy, delta, &mut mdeg, head, invp, perm, qsize, list,
                    marker, maxint, &mut tag,
                );
                // end of -- while ( 1 ) --
                state = 0; // back to top of infinite loop
            }
            3 => {
                // n1000:  line 148
                mmdnum(neqns, perm, invp, qsize);

                /* Adjust from Fortran back to C*/
                // C: xadj++; adjncy++; invp++; perm++; head++; qsize++; list++; marker++;
                // In Rust we never shifted, so nothing to undo.
                return;
            }
            _ => unreachable!(),
        }
    }
}

/// HiGHS: mmd.c::mmdelm (lines 173-292)
#[allow(unused_assignments)] // C locals are initialised before use on every goto path
pub fn mmdelm(
    mdeg_node: idx_t,
    xadj: &mut [idx_t],
    adjncy: &mut [idx_t],
    head: &mut [idx_t],
    forward: &mut [idx_t],
    backward: &mut [idx_t],
    qsize: &mut [idx_t],
    list: &mut [idx_t],
    marker: &mut [idx_t],
    maxint: idx_t,
    tag: idx_t,
) {
    let mut element: idx_t = 0;
    let mut i: idx_t = 0;
    let mut istop: idx_t = 0;
    let mut istart: idx_t = 0;
    let mut j: idx_t = 0;
    let mut jstop: idx_t = 0;
    let mut jstart: idx_t = 0;
    let mut link: idx_t = 0;
    let mut nabor: idx_t = 0;
    let mut node: idx_t = 0;
    let mut npv: idx_t = 0;
    let mut nqnbrs: idx_t = 0;
    let mut nxnode: idx_t = 0;
    let mut pvnode: idx_t = 0;
    let mut rlmt: idx_t = 0;
    let mut rloc: idx_t = 0;
    let mut rnode: idx_t = 0;
    let mut xqnbr: idx_t = 0;

    /* find the reachable set of 'mdeg_node' and */
    /* place it in the data structure.           */
    marker[(mdeg_node - 1) as usize] = tag;
    istart = xadj[(mdeg_node - 1) as usize];
    istop = xadj[(mdeg_node + 1 - 1) as usize] - 1;

    /* 'element' points to the beginning of the list of  */
    /* eliminated nabors of 'mdeg_node', and 'rloc' gives the */
    /* storage location for the next reachable node.   */
    element = 0;
    rloc = istart;
    rlmt = istop;
    i = istart;
    while i <= istop {
        nabor = adjncy[(i - 1) as usize];
        if nabor == 0 {
            break;
        }
        if marker[(nabor - 1) as usize] < tag {
            marker[(nabor - 1) as usize] = tag;
            if forward[(nabor - 1) as usize] < 0 {
                list[(nabor - 1) as usize] = element;
                element = nabor;
            } else {
                adjncy[(rloc - 1) as usize] = nabor;
                rloc += 1;
            }
        }
        i += 1;
    }

    /* merge with reachable nodes from generalized elements. */
    while element > 0 {
        adjncy[(rlmt - 1) as usize] = -element;
        link = element;

        // n400:  line 213
        // Use a labeled loop for the goto n400.
        'n400: loop {
            jstart = xadj[(link - 1) as usize];
            jstop = xadj[(link + 1 - 1) as usize] - 1;
            j = jstart;
            while j <= jstop {
                node = adjncy[(j - 1) as usize];
                link = -node;
                if node < 0 {
                    continue 'n400; // goto n400
                }
                if node == 0 {
                    break;
                }
                if (marker[(node - 1) as usize] < tag) && (forward[(node - 1) as usize] >= 0) {
                    marker[(node - 1) as usize] = tag;
                    /*use storage from eliminated nodes if necessary.*/
                    while rloc >= rlmt {
                        link = -adjncy[(rlmt - 1) as usize];
                        rloc = xadj[(link - 1) as usize];
                        rlmt = xadj[(link + 1 - 1) as usize] - 1;
                    }
                    adjncy[(rloc - 1) as usize] = node;
                    rloc += 1;
                }
                j += 1;
            }
            break;
        }
        element = list[(element - 1) as usize];
    }
    if rloc <= rlmt {
        adjncy[(rloc - 1) as usize] = 0;
    }

    /* for each node in the reachable set, do the following. */
    link = mdeg_node;

    // n1100:  line 239
    // Use a labeled loop for the goto n1100.
    'n1100: loop {
        istart = xadj[(link - 1) as usize];
        istop = xadj[(link + 1 - 1) as usize] - 1;
        i = istart;
        while i <= istop {
            rnode = adjncy[(i - 1) as usize];
            link = -rnode;
            if rnode < 0 {
                continue 'n1100; // goto n1100
            }
            if rnode == 0 {
                return;
            }

            /* 'rnode' is in the degree list structure. */
            pvnode = backward[(rnode - 1) as usize];
            if (pvnode != 0) && (pvnode != (-maxint)) {
                /* then remove 'rnode' from the structure. */
                nxnode = forward[(rnode - 1) as usize];
                if nxnode > 0 {
                    backward[(nxnode - 1) as usize] = pvnode;
                }
                if pvnode > 0 {
                    forward[(pvnode - 1) as usize] = nxnode;
                }
                npv = -pvnode;
                if pvnode < 0 {
                    head[(npv - 1) as usize] = nxnode;
                }
            }

            /* purge inactive quotient nabors of 'rnode'. */
            jstart = xadj[(rnode - 1) as usize];
            jstop = xadj[(rnode + 1 - 1) as usize] - 1;
            xqnbr = jstart;
            j = jstart;
            while j <= jstop {
                nabor = adjncy[(j - 1) as usize];
                if nabor == 0 {
                    break;
                }
                if marker[(nabor - 1) as usize] < tag {
                    adjncy[(xqnbr - 1) as usize] = nabor;
                    xqnbr += 1;
                }
                j += 1;
            }

            /* no active nabor after the purging. */
            nqnbrs = xqnbr - jstart;
            if nqnbrs <= 0 {
                /* merge 'rnode' with 'mdeg_node'. */
                qsize[(mdeg_node - 1) as usize] += qsize[(rnode - 1) as usize];
                qsize[(rnode - 1) as usize] = 0;
                marker[(rnode - 1) as usize] = maxint;
                forward[(rnode - 1) as usize] = -mdeg_node;
                backward[(rnode - 1) as usize] = -maxint;
            } else {
                /* flag 'rnode' for degree update, and  */
                /* add 'mdeg_node' as a nabor of 'rnode'.      */
                forward[(rnode - 1) as usize] = nqnbrs + 1;
                backward[(rnode - 1) as usize] = 0;
                adjncy[(xqnbr - 1) as usize] = mdeg_node;
                xqnbr += 1;
                if xqnbr <= jstop {
                    adjncy[(xqnbr - 1) as usize] = 0;
                }
            }
            i += 1;
        }
        return;
    }
}

/// HiGHS: mmd.c::mmdint (lines 308-332)
#[allow(unused_assignments)] // C locals are initialised before use on every goto path
pub fn mmdint(
    neqns: idx_t,
    xadj: &[idx_t],
    adjncy: &[idx_t],
    head: &mut [idx_t],
    forward: &mut [idx_t],
    backward: &mut [idx_t],
    qsize: &mut [idx_t],
    list: &mut [idx_t],
    marker: &mut [idx_t],
) -> idx_t {
    let mut fnode: idx_t = 0;
    let mut ndeg: idx_t = 0;
    let mut node: idx_t = 0;

    node = 1;
    while node <= neqns {
        head[(node - 1) as usize] = 0;
        qsize[(node - 1) as usize] = 1;
        marker[(node - 1) as usize] = 0;
        list[(node - 1) as usize] = 0;
        node += 1;
    }

    /* initialize the degree doubly linked lists. */
    node = 1;
    while node <= neqns {
        ndeg = xadj[(node + 1 - 1) as usize] - xadj[(node - 1) as usize] + 1;
        fnode = head[(ndeg - 1) as usize];
        forward[(node - 1) as usize] = fnode;
        head[(ndeg - 1) as usize] = node;
        if fnode > 0 {
            backward[(fnode - 1) as usize] = node;
        }
        backward[(node - 1) as usize] = -ndeg;
        node += 1;
    }

    0
}

/// HiGHS: mmd.c::mmdnum (lines 352-396)
#[allow(unused_assignments)] // C locals are initialised before use on every goto path
pub fn mmdnum(neqns: idx_t, perm: &mut [idx_t], invp: &mut [idx_t], qsize: &mut [idx_t]) {
    let mut father: idx_t = 0;
    let mut nextf: idx_t = 0;
    let mut node: idx_t = 0;
    let mut nqsize: idx_t = 0;
    let mut num: idx_t = 0;
    let mut root: idx_t = 0;

    node = 1;
    while node <= neqns {
        nqsize = qsize[(node - 1) as usize];
        if nqsize <= 0 {
            perm[(node - 1) as usize] = invp[(node - 1) as usize];
        }
        if nqsize > 0 {
            perm[(node - 1) as usize] = -invp[(node - 1) as usize];
        }
        node += 1;
    }

    /* for each node which has been merged, do the following. */
    node = 1;
    while node <= neqns {
        if perm[(node - 1) as usize] <= 0 {
            /* trace the merged tree until one which has not */
            /* been merged, call it root.                    */
            father = node;
            while perm[(father - 1) as usize] <= 0 {
                father = -perm[(father - 1) as usize];
            }

            /* number node after root. */
            root = father;
            num = perm[(root - 1) as usize] + 1;
            invp[(node - 1) as usize] = -num;
            perm[(root - 1) as usize] = num;

            /* shorten the merged tree. */
            father = node;
            nextf = -perm[(father - 1) as usize];
            while nextf > 0 {
                perm[(father - 1) as usize] = -root;
                father = nextf;
                nextf = -perm[(father - 1) as usize];
            }
        }
        node += 1;
    }

    /* ready to compute perm. */
    node = 1;
    while node <= neqns {
        num = -invp[(node - 1) as usize];
        invp[(node - 1) as usize] = num;
        perm[(num - 1) as usize] = node;
        node += 1;
    }
}

/// HiGHS: mmd.c::mmdupd (lines 417-598)
#[allow(unused_assignments)] // C locals are initialised before use on every goto path
pub fn mmdupd(
    ehead: idx_t,
    neqns: idx_t,
    xadj: &mut [idx_t],
    adjncy: &mut [idx_t],
    delta: idx_t,
    mdeg: &mut idx_t,
    head: &mut [idx_t],
    forward: &mut [idx_t],
    backward: &mut [idx_t],
    qsize: &mut [idx_t],
    list: &mut [idx_t],
    marker: &mut [idx_t],
    maxint: idx_t,
    tag: &mut idx_t,
) {
    let mut deg: idx_t = 0;
    let mut deg0: idx_t = 0;
    let mut element: idx_t = 0;
    let mut enode: idx_t = 0;
    let mut fnode: idx_t = 0;
    let mut i: idx_t = 0;
    let mut iq2: idx_t = 0;
    let mut istop: idx_t = 0;
    let mut istart: idx_t = 0;
    let mut j: idx_t = 0;
    let mut jstop: idx_t = 0;
    let mut jstart: idx_t = 0;
    let mut link: idx_t = 0;
    let mut mdeg0: idx_t = 0;
    let mut mtag: idx_t = 0;
    let mut nabor: idx_t = 0;
    let mut node: idx_t = 0;
    let mut q2head: idx_t = 0;
    let mut qxhead: idx_t = 0;

    mdeg0 = *mdeg + delta;
    element = ehead;

    // State machine for goto labels:
    // 0 = n100 (top of element loop)
    // 1 = n400 (inner loop over element adjacency, with link-based iteration)
    // 2 = n900 (q2 list processing)
    // 3 = n1000 (inner loop over 2nd element adjacency)
    // 4 = n1500 (switch to qx list)
    // 5 = n1600 (qx list processing)
    // 6 = n1700 (inner loop over eliminated nabor adjacency)
    // 7 = n2100 (update degree)
    // 8 = n2200 (get next enode)
    // 9 = n2300 (get next element)

    let mut state: i32 = 0; // start at n100

    'goto: loop {
        match state {
            0 => {
                // n100:
                if element <= 0 {
                    return;
                }

                /* for each of the newly formed element, do the following. */
                /* reset tag value if necessary.                           */
                mtag = *tag + mdeg0;
                if mtag >= maxint {
                    *tag = 1;
                    i = 1;
                    while i <= neqns {
                        if marker[(i - 1) as usize] < maxint {
                            marker[(i - 1) as usize] = 0;
                        }
                        i += 1;
                    }
                    mtag = *tag + mdeg0;
                }

                /* create two linked lists from nodes associated with 'element': */
                /* one with two nabors (q2head) in the adjacency structure, and the*/
                /* other with more than two nabors (qxhead). also compute 'deg0',*/
                /* number of nodes in this element.                              */
                q2head = 0;
                qxhead = 0;
                deg0 = 0;
                link = element;

                state = 1; // goto n400
            }
            1 => {
                // n400:
                istart = xadj[(link - 1) as usize];
                istop = xadj[(link + 1 - 1) as usize] - 1;
                i = istart;
                while i <= istop {
                    enode = adjncy[(i - 1) as usize];
                    link = -enode;
                    if enode < 0 {
                        state = 1; // goto n400
                        continue 'goto;
                    }
                    if enode == 0 {
                        break;
                    }
                    if qsize[(enode - 1) as usize] != 0 {
                        deg0 += qsize[(enode - 1) as usize];
                        marker[(enode - 1) as usize] = mtag;

                        /*'enode' requires a degree update*/
                        if backward[(enode - 1) as usize] == 0 {
                            /* place either in qxhead or q2head list. */
                            if forward[(enode - 1) as usize] != 2 {
                                list[(enode - 1) as usize] = qxhead;
                                qxhead = enode;
                            } else {
                                list[(enode - 1) as usize] = q2head;
                                q2head = enode;
                            }
                        }
                    }
                    i += 1;
                }

                /* for each node in q2 list, do the following. */
                enode = q2head;
                iq2 = 1;

                state = 2; // goto n900
            }
            2 => {
                // n900:
                if enode <= 0 {
                    state = 4; // goto n1500
                    continue;
                }
                if backward[(enode - 1) as usize] != 0 {
                    state = 8; // goto n2200
                    continue;
                }
                *tag += 1;
                deg = deg0;

                /* identify the other adjacent element nabor. */
                istart = xadj[(enode - 1) as usize];
                nabor = adjncy[(istart - 1) as usize];
                if nabor == element {
                    nabor = adjncy[(istart + 1 - 1) as usize];
                }
                link = nabor;
                if forward[(nabor - 1) as usize] >= 0 {
                    /* nabor is uneliminated, increase degree count. */
                    deg += qsize[(nabor - 1) as usize];
                    state = 7; // goto n2100
                    continue;
                }

                /* the nabor is eliminated. for each node in the 2nd element */
                /* do the following.                                         */
                state = 3; // goto n1000
            }
            3 => {
                // n1000:
                istart = xadj[(link - 1) as usize];
                istop = xadj[(link + 1 - 1) as usize] - 1;
                i = istart;
                while i <= istop {
                    node = adjncy[(i - 1) as usize];
                    link = -node;
                    if node != enode {
                        if node < 0 {
                            state = 3; // goto n1000
                            continue 'goto;
                        }
                        if node == 0 {
                            state = 7; // goto n2100
                            break;
                        }
                        if qsize[(node - 1) as usize] != 0 {
                            if marker[(node - 1) as usize] < *tag {
                                /* 'node' is not yet considered. */
                                marker[(node - 1) as usize] = *tag;
                                deg += qsize[(node - 1) as usize];
                            } else {
                                if backward[(node - 1) as usize] == 0 {
                                    if forward[(node - 1) as usize] == 2 {
                                        /* 'node' is indistinguishable from 'enode'.*/
                                        /* merge them into a new supernode.         */
                                        qsize[(enode - 1) as usize] += qsize[(node - 1) as usize];
                                        qsize[(node - 1) as usize] = 0;
                                        marker[(node - 1) as usize] = maxint;
                                        forward[(node - 1) as usize] = -enode;
                                        backward[(node - 1) as usize] = -maxint;
                                    } else {
                                        /* 'node' is outmatched by 'enode' */
                                        if backward[(node - 1) as usize] == 0 {
                                            backward[(node - 1) as usize] = -maxint;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    i += 1;
                }
                state = 7; // goto n2100
            }
            4 => {
                // n1500:
                /* for each 'enode' in the 'qx' list, do the following. */
                enode = qxhead;
                iq2 = 0;

                state = 5; // goto n1600
            }
            5 => {
                // n1600:
                if enode <= 0 {
                    state = 9; // goto n2300
                    continue;
                }
                if backward[(enode - 1) as usize] != 0 {
                    state = 8; // goto n2200
                    continue;
                }
                *tag += 1;
                deg = deg0;

                /*for each unmarked nabor of 'enode', do the following.*/
                istart = xadj[(enode - 1) as usize];
                istop = xadj[(enode + 1 - 1) as usize] - 1;
                i = istart;
                while i <= istop {
                    nabor = adjncy[(i - 1) as usize];
                    if nabor == 0 {
                        break;
                    }
                    if marker[(nabor - 1) as usize] < *tag {
                        marker[(nabor - 1) as usize] = *tag;
                        link = nabor;
                        if forward[(nabor - 1) as usize] >= 0 {
                            /*if uneliminated, include it in deg count.*/
                            deg += qsize[(nabor - 1) as usize];
                        } else {
                            // n1700:
                            // The goto n1700 from within the inner loop means
                            // we restart the inner loop with a new link.
                            // Use a labeled loop to faithfully translate the goto.
                            'n1700: loop {
                                jstart = xadj[(link - 1) as usize];
                                jstop = xadj[(link + 1 - 1) as usize] - 1;
                                j = jstart;
                                while j <= jstop {
                                    node = adjncy[(j - 1) as usize];
                                    link = -node;
                                    if node < 0 {
                                        continue 'n1700;
                                    }
                                    if node == 0 {
                                        break;
                                    }
                                    if marker[(node - 1) as usize] < *tag {
                                        marker[(node - 1) as usize] = *tag;
                                        deg += qsize[(node - 1) as usize];
                                    }
                                    j += 1;
                                }
                                break;
                            }
                        }
                    }
                    i += 1;
                }

                state = 7; // goto n2100
            }
            7 => {
                // n2100:
                /* update external degree of 'enode' in degree structure, */
                /* and '*mdeg' if necessary.                     */
                deg = deg - qsize[(enode - 1) as usize] + 1;
                fnode = head[(deg - 1) as usize];
                forward[(enode - 1) as usize] = fnode;
                backward[(enode - 1) as usize] = -deg;
                if fnode > 0 {
                    backward[(fnode - 1) as usize] = enode;
                }
                head[(deg - 1) as usize] = enode;
                if deg < *mdeg {
                    *mdeg = deg;
                }

                state = 8; // goto n2200
            }
            8 => {
                // n2200:
                /* get next enode in current element. */
                enode = list[(enode - 1) as usize];
                if iq2 == 1 {
                    state = 2; // goto n900
                } else {
                    state = 5; // goto n1600
                }
            }
            9 => {
                // n2300:
                /* get next element in the list. */
                *tag = mtag;
                element = list[(element - 1) as usize];
                state = 0; // goto n100
            }
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn genmmd_path_graph() {
        // 1-based CSR of the path 1-2-3 (as MMDOrder passes it)
        let mut xadj = vec![1, 2, 4, 5];
        let mut adjncy = vec![2, 1, 3, 2];
        let n = 3;
        let (mut invp, mut perm) = (vec![0; n + 5], vec![0; n + 5]);
        let (mut head, mut qsize, mut list, mut marker) = (
            vec![0; n + 5],
            vec![0; n + 5],
            vec![0; n + 5],
            vec![0; n + 5],
        );
        let mut nofsub = 0;
        genmmd(
            n as idx_t,
            &mut xadj,
            &mut adjncy,
            &mut invp,
            &mut perm,
            1,
            &mut head,
            &mut qsize,
            &mut list,
            &mut marker,
            crate::ipm::hipo::orderings::metis::IDX_MAX,
            &mut nofsub,
        );
        let mut p: Vec<idx_t> = perm[..n].to_vec();
        p.sort();
        assert_eq!(p, vec![1, 2, 3]);
    }

    /// Hand-built state reaching every `mmdupd` label: n400 and n1000 link
    /// chains (negative entries), a q2 merge, an outmatched node, n1700 link
    /// chains, a second element, and (with `tag0` near `maxint`) the tag reset.
    /// Expected values were produced by compiling mmd.c `mmdupd` unmodified.
    fn run_mmdupd(tag0: idx_t, maxint: idx_t) -> (idx_t, idx_t, [Vec<idx_t>; 6]) {
        let mut xadj = vec![1, 3, 5, 7, 9, 11, 14, 17, 19, 21, 23];
        let mut adjncy = vec![
            2, -3, 1, 9, 4, -6, 1, 9, 6, 0, 5, 7, 0, 1, 8, 0, 7, 0, 2, -10, 4, 0,
        ];
        let mut forward = vec![-1, 2, -2, 2, 3, -3, 5, 0, -4, -5];
        let mut backward = vec![0; 10];
        let mut qsize = vec![1; 10];
        let mut list = vec![0; 10];
        let mut marker = vec![0; 10];
        let mut head = vec![0; 10];
        list[0] = 9;
        marker[7] = maxint - 5;
        marker[2] = maxint;
        let (mut mdeg, mut tag) = (3, tag0);
        mmdupd(
            1,
            10,
            &mut xadj,
            &mut adjncy,
            1,
            &mut mdeg,
            &mut head,
            &mut forward,
            &mut backward,
            &mut qsize,
            &mut list,
            &mut marker,
            maxint,
            &mut tag,
        );
        (mdeg, tag, [forward, backward, qsize, list, marker, head])
    }

    #[test]
    fn mmdupd_all_labels() {
        let (mdeg, tag, [forward, backward, qsize, list, marker, head]) = run_mmdupd(1, 1000);
        assert_eq!((mdeg, tag), (3, 9));
        assert_eq!(forward, vec![-1, -4, -2, 0, 7, -3, 0, 0, -4, -5]);
        assert_eq!(backward, vec![0, -1000, 0, -3, -4, 0, 5, 0, 0, 0]);
        assert_eq!(qsize, vec![1, 0, 1, 2, 1, 1, 1, 1, 1, 1]);
        assert_eq!(list, vec![9, 0, 0, 2, 0, 0, 5, 0, 0, 0]);
        assert_eq!(marker, vec![3, 1000, 1000, 9, 5, 4, 5, 995, 0, 0]);
        assert_eq!(head, vec![0, 0, 4, 5, 0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn mmdupd_tag_reset() {
        let (mdeg, tag, [forward, backward, qsize, list, marker, head]) = run_mmdupd(997, 1000);
        assert_eq!((mdeg, tag), (3, 9));
        assert_eq!(forward, vec![-1, -4, -2, 0, 0, -3, 0, 0, -4, -5]);
        assert_eq!(backward, vec![0, -1000, 0, -3, -4, 0, -5, 0, 0, 0]);
        assert_eq!(qsize, vec![1, 0, 1, 2, 1, 1, 1, 1, 1, 1]);
        assert_eq!(list, vec![9, 0, 0, 2, 0, 0, 5, 0, 0, 0]);
        assert_eq!(marker, vec![3, 1000, 1000, 9, 5, 4, 5, 3, 0, 0]);
        assert_eq!(head, vec![0, 0, 4, 5, 7, 0, 0, 0, 0, 0]);
    }

    /// Node 2 merged into 1, node 4 merged into 2 (a chain the shortening
    /// loop rewrites). Expected values from mmd.c `mmdnum` compiled unmodified.
    #[test]
    fn mmdnum_merged_chain() {
        let mut qsize = vec![3, 0, 1, 0, 1];
        let mut invp = vec![-1, -1, -4, -2, -5];
        let mut perm = vec![0; 5];
        mmdnum(5, &mut perm, &mut invp, &mut qsize);
        assert_eq!(perm, vec![1, 2, 4, 3, 5]);
        assert_eq!(invp, vec![1, 2, 4, 3, 5]);
    }
    /// Runs `genmmd` with every work vector sized `n + 5`, as `MMDOrder` allocates them.
    /// Returns `ncsub`, then invp, perm, head, qsize, list, marker, and the destroyed adjncy.
    fn run_genmmd(
        mut xadj: Vec<idx_t>,
        mut adjncy: Vec<idx_t>,
        delta: idx_t,
        maxint: idx_t,
    ) -> (idx_t, [Vec<idx_t>; 7]) {
        let n = xadj.len() - 1;
        let mut v: [Vec<idx_t>; 6] = std::array::from_fn(|_| vec![0; n + 5]);
        let [invp, perm, head, qsize, list, marker] = &mut v;
        let mut ncsub = 0;
        genmmd(
            n as idx_t,
            &mut xadj,
            &mut adjncy,
            invp,
            perm,
            delta,
            head,
            qsize,
            list,
            marker,
            maxint,
            &mut ncsub,
        );
        let [invp, perm, head, qsize, list, marker] = v;
        (ncsub, [invp, perm, head, qsize, list, marker, adjncy])
    }

    /// Random 17-node graph with an isolated node (10), delta = 1 and maxint = 46: reaches
    /// the genmmd tag reset (mmd.c 128-133), mmdelm's n400 element chain (215), storage
    /// reuse (223-227), n1100 link chain (243) and the qsize merge (268-273).
    /// Expected values from mmd.c `genmmd` compiled unmodified (gcc, ASan/UBSan clean).
    #[test]
    fn genmmd_random_tag_reset() {
        let xadj = vec![
            1, 5, 8, 10, 15, 17, 21, 22, 26, 30, 30, 31, 33, 34, 36, 38, 41, 45,
        ];
        let adjncy = vec![
            3, 5, 6, 9, 4, 6, 7, 1, 8, 2, 6, 14, 16, 17, 1, 9, 1, 2, 4, 14, 2, 3, 11, 16, 17, 1, 5,
            15, 16, 8, 15, 17, 17, 4, 6, 9, 12, 4, 8, 9, 4, 8, 12, 13,
        ];
        let (ncsub, [invp, perm, head, qsize, list, marker, adjncy]) =
            run_genmmd(xadj, adjncy, 1, 46);
        assert_eq!(ncsub, 29);
        let z5 = [0; 5];
        let ext = |a: &[idx_t]| [a, &z5[..]].concat();
        assert_eq!(
            invp,
            ext(&[11, 9, 8, 16, 7, 14, 4, 15, 17, 1, 3, 10, 2, 6, 5, 12, 13])
        );
        assert_eq!(
            perm,
            ext(&[10, 13, 11, 7, 15, 14, 5, 3, 2, 12, 1, 16, 17, 6, 8, 4, 9])
        );
        assert_eq!(head, vec![0; 22]);
        assert_eq!(
            qsize,
            ext(&[1, 1, 1, 0, 1, 1, 1, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1])
        );
        assert_eq!(
            list,
            ext(&[2, 14, 5, 9, 0, 17, 0, 0, 8, 0, 0, 13, 0, 0, 0, 11, 16])
        );
        assert_eq!(
            marker,
            ext(&[
                28, 28, 1, 46, 1, 43, 0, 43, 46, 46, 43, 27, 27, 28, 0, 43, 43
            ])
        );
        assert_eq!(
            adjncy,
            vec![
                6, 8, 9, 0, 4, 6, 0, 1, 8, 16, 17, 6, 0, 0, 1, 9, 8, 9, 4, 0, 2, 9, 4, 0, -11, 16,
                17, 6, 0, 8, 17, -15, 17, 4, 6, 9, 0, 4, 8, 9, 4, 8, 9, 0,
            ]
        );
    }

    /// delta = -1 (no multiple elimination: mmdupd after every node, mmd.c 136-137) on a
    /// 7-node graph that reaches the n400 element chain and two qsize merges.
    /// Expected values from mmd.c `genmmd` compiled unmodified.
    #[test]
    fn genmmd_negative_delta() {
        let xadj = vec![1, 5, 6, 8, 10, 12, 12, 15];
        let adjncy = vec![2, 3, 4, 7, 1, 1, 5, 1, 7, 3, 7, 1, 4, 5];
        let maxint = crate::ipm::hipo::orderings::metis::IDX_MAX;
        let (ncsub, [invp, perm, head, qsize, list, marker, adjncy]) =
            run_genmmd(xadj, adjncy, -1, maxint);
        assert_eq!(ncsub, 7);
        let z5 = [0; 5];
        let ext = |a: &[idx_t]| [a, &z5[..]].concat();
        assert_eq!(invp, ext(&[5, 2, 4, 6, 3, 1, 7]));
        assert_eq!(perm, ext(&[6, 2, 5, 3, 1, 4, 7]));
        assert_eq!(head, vec![0; 12]);
        assert_eq!(qsize, ext(&[3, 1, 1, 0, 1, 1, 0]));
        assert_eq!(list, ext(&[0, 0, 2, 0, 0, 0, 0]));
        assert_eq!(marker, ext(&[10, 10, 10, maxint, 7, maxint, maxint]));
        assert_eq!(adjncy, vec![4, 7, 0, -2, 1, 1, -5, 1, 7, 7, 0, 4, 3, 0]);
    }
}
