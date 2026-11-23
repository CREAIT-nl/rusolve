//! HiGHS source: extern/amd/amd_postorder.c, extern/amd/amd_post_tree.c
//!
//! Translated from SuiteSparse AMD as vendored by HiGHS in `extern/amd` (BSD-3-Clause,
//! Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff).
//! See THIRD_PARTY_NOTICES.md.
//!
//! Postorder of the assembly tree (children ordered so that the child with the
//! largest front is visited last) and the non-recursive depth-first numbering.

use super::EMPTY;
use super::amd_int;

/// HiGHS: amd_postorder.c::Highs_amd_postorder (lines 15-140)
pub fn Highs_amd_postorder(
    nn: amd_int,
    Parent: &[amd_int],
    Nv: &[amd_int],
    Fsize: &[amd_int],
    Order: &mut [amd_int],
    Child: &mut [amd_int],
    Sibling: &mut [amd_int],
    Stack: &mut [amd_int],
) {
    let mut i: amd_int;
    let mut j: amd_int;
    let mut k: amd_int;
    let mut parent: amd_int;
    let mut frsize: amd_int;
    let mut f: amd_int;
    let mut fprev: amd_int;
    let mut maxfrsize: amd_int;
    let mut bigfprev: amd_int;
    let mut bigf: amd_int;
    let mut fnext: amd_int;

    j = 0;
    while j < nn {
        Child[j as usize] = EMPTY;
        Sibling[j as usize] = EMPTY;
        j += 1;
    }

    /* --------------------------------------------------------------------- */
    /* place the children in link lists - bigger elements tend to be last */
    /* --------------------------------------------------------------------- */

    j = nn - 1;
    while j >= 0 {
        if Nv[j as usize] > 0 {
            /* this is an element */
            parent = Parent[j as usize];
            if parent != EMPTY {
                /* place the element in link list of the children its parent */
                /* bigger elements will tend to be at the end of the list */
                Sibling[j as usize] = Child[parent as usize];
                Child[parent as usize] = j;
            }
        }
        j -= 1;
    }

    /* --------------------------------------------------------------------- */
    /* place the largest child last in the list of children for each node */
    /* --------------------------------------------------------------------- */

    i = 0;
    while i < nn {
        if Nv[i as usize] > 0 && Child[i as usize] != EMPTY {
            /* find the biggest element in the child list */
            fprev = EMPTY;
            maxfrsize = EMPTY;
            bigfprev = EMPTY;
            bigf = EMPTY;
            f = Child[i as usize];
            while f != EMPTY {
                frsize = Fsize[f as usize];
                if frsize >= maxfrsize {
                    /* this is the biggest seen so far */
                    maxfrsize = frsize;
                    bigfprev = fprev;
                    bigf = f;
                }
                fprev = f;
                f = Sibling[f as usize];
            }

            fnext = Sibling[bigf as usize];

            if fnext != EMPTY {
                /* if fnext is EMPTY then bigf is already at the end of list */

                if bigfprev == EMPTY {
                    /* delete bigf from the element of the list */
                    Child[i as usize] = fnext;
                } else {
                    /* delete bigf from the middle of the list */
                    Sibling[bigfprev as usize] = fnext;
                }

                /* put bigf at the end of the list */
                Sibling[bigf as usize] = EMPTY;

                Sibling[fprev as usize] = bigf;
            }
        }
        i += 1;
    }

    /* --------------------------------------------------------------------- */
    /* postorder the assembly tree */
    /* --------------------------------------------------------------------- */

    i = 0;
    while i < nn {
        Order[i as usize] = EMPTY;
        i += 1;
    }

    k = 0;

    i = 0;
    while i < nn {
        if Parent[i as usize] == EMPTY && Nv[i as usize] > 0 {
            k = Highs_amd_post_tree(i, k, Child, Sibling, Order, Stack);
        }
        i += 1;
    }
}

/// HiGHS: amd_post_tree.c::Highs_amd_post_tree (lines 15-105)
pub fn Highs_amd_post_tree(
    root: amd_int,
    k: amd_int,
    Child: &mut [amd_int],
    Sibling: &[amd_int],
    Order: &mut [amd_int],
    Stack: &mut [amd_int],
) -> amd_int {
    let mut f: amd_int;
    let mut head: amd_int;
    let mut h: amd_int;
    let mut i: amd_int;
    let mut k = k;

    // HiGHS: amd_post_tree.c lines 36-49: #if 0 recursive version dropped

    /* --------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* --------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0;
    Stack[0] = root;

    while head >= 0 {
        /* get head of stack */
        i = Stack[head as usize];

        if Child[i as usize] != EMPTY {
            /* the children of i are not yet ordered */
            /* push each child onto the stack in reverse order */
            /* so that small ones at the head of the list get popped first */
            /* and the biggest one at the end of the list gets popped last */
            f = Child[i as usize];
            while f != EMPTY {
                head += 1;
                f = Sibling[f as usize];
            }
            h = head;

            f = Child[i as usize];
            while f != EMPTY {
                Stack[h as usize] = f;
                h -= 1;
                f = Sibling[f as usize];
            }

            /* delete child list so that i gets ordered next time we see it */
            Child[i as usize] = EMPTY;
        } else {
            /* the children of i (if there were any) are already ordered */
            /* remove i from the stack and order it.  Front i is kth front */
            head -= 1;

            Order[i as usize] = k;
            k += 1;
        }
    }
    k
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn postorder_chain() {
        // chain 0 -> 1 -> 2 (2 is the root)
        let parent = [1, 2, crate::ipm::hipo::orderings::amd::EMPTY];
        let nv = [1, 1, 1];
        let fsize = [1, 1, 1];
        let mut order = [0; 3];
        let (mut child, mut sibling, mut stack) = ([0; 3], [0; 3], [0; 3]);
        Highs_amd_postorder(
            3,
            &parent,
            &nv,
            &fsize,
            &mut order,
            &mut child,
            &mut sibling,
            &mut stack,
        );
        assert_eq!(order, [0, 1, 2]);
    }

    fn run(parent: &[amd_int], nv: &[amd_int], fsize: &[amd_int]) -> Vec<amd_int> {
        let nn = parent.len();
        let mut order = vec![0; nn];
        let (mut child, mut sibling, mut stack) = (vec![0; nn], vec![0; nn], vec![0; nn]);
        Highs_amd_postorder(
            nn as amd_int,
            parent,
            nv,
            fsize,
            &mut order,
            &mut child,
            &mut sibling,
            &mut stack,
        );
        order
    }

    #[test]
    fn postorder_biggest_child_at_head_moved_last() {
        // children 0,1,2 of root 3; child 0 (list head) is biggest:
        // amd_postorder.c bigfprev == EMPTY branch, list becomes 1,2,0
        let order = run(&[3, 3, 3, EMPTY], &[1, 1, 1, 1], &[5, 1, 1, 1]);
        assert_eq!(order, [2, 0, 1, 3]);
    }

    #[test]
    fn postorder_biggest_child_in_middle_moved_last() {
        // child 1 (middle of list 0,1,2) is biggest: Sibling[bigfprev] = fnext
        // branch, list becomes 0,2,1
        let order = run(&[3, 3, 3, EMPTY], &[1, 1, 1, 1], &[1, 5, 1, 1]);
        assert_eq!(order, [0, 2, 1, 3]);
    }

    #[test]
    fn postorder_multiple_roots_and_non_nodes() {
        // roots 0 and 1; node 2 child of 0; node 3 has Nv == 0 (not a node),
        // so it is neither linked nor ordered and keeps Order[3] == EMPTY
        let order = run(&[EMPTY, EMPTY, 0, 0], &[1, 1, 1, 0], &[1, 1, 1, 1]);
        assert_eq!(order, [1, 2, 0, EMPTY]);
    }

    #[test]
    fn post_tree_returns_next_k() {
        // single root 0 with children 1,2 (list 1,2), numbering starts at k = 4
        let mut child = [1, EMPTY, EMPTY];
        let sibling = [EMPTY, 2, EMPTY];
        let mut order = [EMPTY; 3];
        let mut stack = [0; 3];
        let k = Highs_amd_post_tree(0, 4, &mut child, &sibling, &mut order, &mut stack);
        assert_eq!(k, 7);
        assert_eq!(order, [6, 4, 5]);
        assert_eq!(child, [EMPTY, EMPTY, EMPTY]);
    }
}
