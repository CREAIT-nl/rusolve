//! HiGHS source: extern/metis/GKlib/gk_mkpqueue.h (GK_MKPQUEUE(rpq, rpq_t, rkv_t, real_t, idx_t, rkvmalloc, REAL_MAX, key_gt), libmetis/gklib.c line 32)
//!
//! Translated from METIS 5.2.1 as vendored by HiGHS in `extern/metis` (Apache License 2.0,
//! Copyright 1995-2013, Regents of the University of Minnesota). See THIRD_PARTY_NOTICES.md.
//!
//! The max-heap priority queue of `real_t` keys (`key_gt` = `>`), with a
//! `locator` array for update/delete by node. Heap sift order is ported
//! exactly: ties between equal keys are resolved by the C comparisons, which
//! decide the FM move order.

use super::{REAL_MAX, idx_t, real_t, rkv_t};

/// HiGHS: gk_struct.h `GK_MKPQUEUE_T(rpq_t, rkv_t)`
#[derive(Clone, Debug, Default)]
pub struct rpq_t {
    pub nnodes: usize,
    pub maxnodes: usize,
    pub heap: Vec<rkv_t>,
    pub locator: Vec<isize>,
}

/// HiGHS: gk_mkpqueue.h::rpqCreate (lines 19-27)
pub fn rpqCreate(maxnodes: usize) -> rpq_t {
    let mut queue = rpq_t::default();
    rpqInit(&mut queue, maxnodes);
    queue
}

/// HiGHS: gk_mkpqueue.h::rpqInit (lines 33-42)
/// `locator` is filled with -1 (`gk_idxsmalloc(maxnodes, -1)`).
pub fn rpqInit(queue: &mut rpq_t, maxnodes: usize) {
    queue.nnodes = 0;
    queue.maxnodes = maxnodes;

    queue.heap = vec![rkv_t::default(); maxnodes];
    queue.locator = vec![-1isize; maxnodes];
}

/// HiGHS: gk_mkpqueue.h::rpqReset (lines 48-57)
pub fn rpqReset(queue: &mut rpq_t) {
    let locator = &mut queue.locator;
    let heap = &queue.heap;

    let mut i: isize = queue.nnodes as isize - 1;
    while i >= 0 {
        locator[heap[i as usize].val as usize] = -1;
        i -= 1;
    }
    queue.nnodes = 0;
}

/// HiGHS: gk_mkpqueue.h::rpqFree (lines 63-69)
pub fn rpqFree(queue: &mut rpq_t) {
    // C: if (queue == NULL) return; — Rust ownership makes this impossible.
    // gk_free(&heap) / gk_free(&locator): release the storage (C sets them to NULL).
    queue.heap = Vec::new();
    queue.locator = Vec::new();
    queue.maxnodes = 0;
}

/// HiGHS: gk_mkpqueue.h::rpqDestroy (lines 76-81)
pub fn rpqDestroy(queue: rpq_t) {
    // C: if (queue == NULL) return; — Rust ownership makes this impossible.
    // rpqFree(queue); gk_free(&queue); — Drop handles deallocation.
    drop(queue);
}

/// HiGHS: gk_mkpqueue.h::rpqLength (lines 87-90)
pub fn rpqLength(queue: &rpq_t) -> usize {
    queue.nnodes
}

/// HiGHS: gk_mkpqueue.h::rpqInsert (lines 96-118)
pub fn rpqInsert(queue: &mut rpq_t, node: idx_t, key: real_t) -> i32 {
    let mut i: isize;
    let mut j: isize;
    let locator = &mut queue.locator;
    let heap = &mut queue.heap;

    i = queue.nnodes as isize;
    queue.nnodes += 1;
    while i > 0 {
        j = (i - 1) >> 1;
        if key > heap[j as usize].key {
            heap[i as usize] = heap[j as usize];
            locator[heap[i as usize].val as usize] = i;
            i = j;
        } else {
            break;
        }
    }
    heap[i as usize].key = key;
    heap[i as usize].val = node;
    locator[node as usize] = i;

    0
}

/// HiGHS: gk_mkpqueue.h::rpqDelete (lines 124-179)
pub fn rpqDelete(queue: &mut rpq_t, node: idx_t) -> i32 {
    let mut i: isize;
    let mut j: isize;
    let mut nnodes: usize;
    let mut newkey: real_t;
    let mut oldkey: real_t;
    let locator = &mut queue.locator;
    let heap = &mut queue.heap;

    i = locator[node as usize];
    locator[node as usize] = -1;

    queue.nnodes -= 1;
    if queue.nnodes > 0 && heap[queue.nnodes].val != node {
        let mut node = heap[queue.nnodes].val;
        newkey = heap[queue.nnodes].key;
        oldkey = heap[i as usize].key;

        if newkey > oldkey {
            // Filter-up
            while i > 0 {
                j = (i - 1) >> 1;
                if newkey > heap[j as usize].key {
                    heap[i as usize] = heap[j as usize];
                    locator[heap[i as usize].val as usize] = i;
                    i = j;
                } else {
                    break;
                }
            }
        } else {
            // Filter down
            nnodes = queue.nnodes;
            loop {
                j = (i << 1) + 1;
                if j >= nnodes as isize {
                    break;
                }
                if heap[j as usize].key > newkey {
                    if (j + 1) < nnodes as isize
                        && heap[(j + 1) as usize].key > heap[j as usize].key
                    {
                        j += 1;
                    }
                    heap[i as usize] = heap[j as usize];
                    locator[heap[i as usize].val as usize] = i;
                    i = j;
                } else if (j + 1) < nnodes as isize && heap[(j + 1) as usize].key > newkey {
                    j += 1;
                    heap[i as usize] = heap[j as usize];
                    locator[heap[i as usize].val as usize] = i;
                    i = j;
                } else {
                    break;
                }
            }
        }

        heap[i as usize].key = newkey;
        heap[i as usize].val = node;
        locator[node as usize] = i;
    }

    0
}

/// HiGHS: gk_mkpqueue.h::rpqUpdate (lines 185-236)
pub fn rpqUpdate(queue: &mut rpq_t, node: idx_t, newkey: real_t) {
    let mut i: isize;
    let mut j: isize;
    let mut nnodes: usize;

    let locator = &mut queue.locator;
    let heap = &mut queue.heap;

    let oldkey: real_t = heap[locator[node as usize] as usize].key;
    if !(newkey > oldkey) && !(oldkey > newkey) {
        return;
    }

    i = locator[node as usize];

    if newkey > oldkey {
        // Filter-up
        while i > 0 {
            j = (i - 1) >> 1;
            if newkey > heap[j as usize].key {
                heap[i as usize] = heap[j as usize];
                locator[heap[i as usize].val as usize] = i;
                i = j;
            } else {
                break;
            }
        }
    } else {
        // Filter down
        nnodes = queue.nnodes;
        loop {
            j = (i << 1) + 1;
            if j >= nnodes as isize {
                break;
            }
            if heap[j as usize].key > newkey {
                if (j + 1) < nnodes as isize && heap[(j + 1) as usize].key > heap[j as usize].key {
                    j += 1;
                }
                heap[i as usize] = heap[j as usize];
                locator[heap[i as usize].val as usize] = i;
                i = j;
            } else if (j + 1) < nnodes as isize && heap[(j + 1) as usize].key > newkey {
                j += 1;
                heap[i as usize] = heap[j as usize];
                locator[heap[i as usize].val as usize] = i;
                i = j;
            } else {
                break;
            }
        }
    }

    heap[i as usize].key = newkey;
    heap[i as usize].val = node;
    locator[node as usize] = i;
}

/// HiGHS: gk_mkpqueue.h::rpqGetTop (lines 243-290)
/// Returns -1 when empty.
pub fn rpqGetTop(queue: &mut rpq_t) -> idx_t {
    let mut i: isize;
    let mut j: isize;
    let locator = &mut queue.locator;
    let heap = &mut queue.heap;

    let mut node: idx_t;
    let mut key: real_t;

    if queue.nnodes == 0 {
        return -1;
    }

    queue.nnodes -= 1;

    let vtx: idx_t = heap[0].val;
    locator[vtx as usize] = -1;

    i = queue.nnodes as isize;
    if i > 0 {
        key = heap[i as usize].key;
        node = heap[i as usize].val;
        i = 0;
        loop {
            j = 2 * i + 1;
            if j >= queue.nnodes as isize {
                break;
            }
            if heap[j as usize].key > key {
                if (j + 1) < queue.nnodes as isize
                    && heap[(j + 1) as usize].key > heap[j as usize].key
                {
                    j = j + 1;
                }
                heap[i as usize] = heap[j as usize];
                locator[heap[i as usize].val as usize] = i;
                i = j;
            } else if (j + 1) < queue.nnodes as isize && heap[(j + 1) as usize].key > key {
                j = j + 1;
                heap[i as usize] = heap[j as usize];
                locator[heap[i as usize].val as usize] = i;
                i = j;
            } else {
                break;
            }
        }

        heap[i as usize].key = key;
        heap[i as usize].val = node;
        locator[node as usize] = i;
    }

    vtx
}

/// HiGHS: gk_mkpqueue.h::rpqSeeTopVal (lines 297-300)
pub fn rpqSeeTopVal(queue: &rpq_t) -> idx_t {
    if queue.nnodes == 0 {
        -1
    } else {
        queue.heap[0].val
    }
}

/// HiGHS: gk_mkpqueue.h::rpqSeeTopKey (lines 307-310)
/// Returns `REAL_MAX` when empty.
pub fn rpqSeeTopKey(queue: &rpq_t) -> real_t {
    if queue.nnodes == 0 {
        REAL_MAX
    } else {
        queue.heap[0].key
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn insert_gettop_order() {
        let mut q = rpqCreate(8);
        rpqInsert(&mut q, 3, 1.0);
        rpqInsert(&mut q, 5, 4.0);
        rpqInsert(&mut q, 1, 2.0);
        rpqUpdate(&mut q, 3, 9.0);
        assert_eq!(rpqSeeTopKey(&q), 9.0);
        assert_eq!(rpqGetTop(&mut q), 3);
        rpqDelete(&mut q, 1);
        assert_eq!(rpqGetTop(&mut q), 5);
        assert_eq!(rpqGetTop(&mut q), -1);
    }
}
