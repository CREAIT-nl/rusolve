//! HiGHS source: ipm/ipx/sparse_utils.h, sparse_utils.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;

/// HiGHS: ipx/sparse_utils.cc::DepthFirstSearch (lines 6-38)
///
/// `colmap` may be `nullptr` (identity).
#[allow(clippy::too_many_arguments)]
pub fn depth_first_search(
    istart: Int,
    Ap: &[Int],
    Ai: &[Int],
    colmap: Option<&[Int]>,
    top: Int,
    istack: &mut [Int],
    marked: &mut [Int],
    marker: Int,
    work: &mut [Int],
) -> Int {
    debug_assert!(marked[istart as usize] != marker);
    let pstack = work;
    let mut top = top;
    let mut head: Int = 0;
    istack[head as usize] = istart;
    while head >= 0 {
        let i = istack[head as usize];
        let j = match colmap {
            Some(c) => c[i as usize],
            None => i,
        };
        if marked[i as usize] != marker {
            marked[i as usize] = marker;
            pstack[head as usize] = if j >= 0 { Ap[j as usize] } else { 0 };
        }
        let mut done = true;
        let pbeg = pstack[head as usize];
        let pend = if j >= 0 { Ap[(j + 1) as usize] } else { 0 };
        let mut p = pbeg;
        while p < pend {
            let inext = Ai[p as usize];
            if marked[inext as usize] == marker {
                p += 1;
                continue;
            }
            pstack[head as usize] = p + 1;
            head += 1;
            istack[head as usize] = inext;
            done = false;
            break;
        }
        if done {
            head -= 1;
            top -= 1;
            istack[top as usize] = i;
        }
    }
    top
}

/// HiGHS: ipx/sparse_utils.cc::AugmentingPath (lines 40-90)
#[allow(clippy::too_many_arguments)]
pub fn augmenting_path(
    jstart: Int,
    Ap: &[Int],
    Ai: &[Int],
    jmatch: &mut [Int],
    cheap: &mut [Int],
    marked: &mut [Int],
    work: &mut [Int],
    work2: &mut [Int],
    work3: &mut [Int],
) -> bool {
    // istack holds row indices without duplicates.
    // jstack holds jstart and up to m indices from jmatch.
    // pstack holds the corresponding column pointers.
    let istack = work;
    let jstack = work2;
    let pstack = work3;
    let mut found = false;
    let mut head: Int = 0;
    jstack[head as usize] = jstart;
    while head >= 0 {
        let j = jstack[head as usize];
        if marked[j as usize] != jstart {
            // first time j visited in this path
            marked[j as usize] = jstart;
            let mut i: Int = 0;
            let mut p = cheap[j as usize];
            while p < Ap[(j + 1) as usize] && !found {
                i = Ai[p as usize];
                found = jmatch[i as usize] == -1;
                p += 1;
            }
            cheap[j as usize] = p;
            if found {
                istack[head as usize] = i;
                break;
            }
            pstack[head as usize] = Ap[j as usize]; // start depth-first search from j
        }
        let pbeg = pstack[head as usize];
        let pend = Ap[(j + 1) as usize];
        let mut p = pbeg;
        while p < pend {
            let i = Ai[p as usize];
            if jmatch[i as usize] < -1 {
                // row to be ignored
                p += 1;
                continue;
            }
            debug_assert!(jmatch[i as usize] >= 0);
            if marked[jmatch[i as usize] as usize] == jstart {
                p += 1;
                continue;
            }
            pstack[head as usize] = p + 1;
            istack[head as usize] = i;
            head += 1;
            jstack[head as usize] = jmatch[i as usize]; // continue augmenting path here
            break;
        }
        if p == Ap[(j + 1) as usize] {
            head -= 1;
        }
    }
    if found {
        let mut p = head;
        while p >= 0 {
            jmatch[istack[p as usize] as usize] = jstack[p as usize];
            p -= 1;
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dfs_postorder_on_chain() {
        // Graph 0 -> 1 -> 2 in CSC adjacency form: column 0 holds row 1, column 1 holds row 2.
        let ap = [0, 1, 2, 2];
        let ai = [1, 2];
        let mut istack = [0; 3];
        let mut marked = [-1; 3];
        let mut work = [0; 3];
        let top = depth_first_search(0, &ap, &ai, None, 3, &mut istack, &mut marked, 7, &mut work);
        assert_eq!(top, 0);
        assert_eq!(istack, [0, 1, 2]);
    }

    #[test]
    fn augmenting_path_direct_match() {
        // Column 0 has row 0, which is free.
        let ap = [0, 1, 2];
        let ai = [0, 1];
        let mut jmatch = [-1, -1];
        let mut cheap = [0, 1];
        let mut marked = [-1, -1];
        let (mut w1, mut w2, mut w3) = ([0; 2], [0; 2], [0; 2]);
        let found = augmenting_path(
            0,
            &ap,
            &ai,
            &mut jmatch,
            &mut cheap,
            &mut marked,
            &mut w1,
            &mut w2,
            &mut w3,
        );
        assert!(found);
        assert_eq!(jmatch, [0, -1]);
    }
}
