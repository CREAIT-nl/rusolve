//! HiGHS source: extern/rcm/rcm.cpp (lines 1-12, 347-566, 758-879), extern/rcm/rcm.h
//!
//! Translated from RCM as vendored by HiGHS in `extern/rcm` (MIT License, Copyright (c)
//! 2005 Alan George, Joseph Liu, John Burkardt). See THIRD_PARTY_NOTICES.md.
//!
//! Reverse Cuthill-McKee (George & Liu, Burkardt's C++). The C++ keeps the
//! Fortran 1-based node numbers (`offset = 1`): node `i` lives at `a[i - 1]`
//! and `adj_row` / `adj` are read with `+ offset`. The Rust indexes the same
//! way through `(x - 1) as usize`. The C++ passes `perm + num - 1` (an
//! interior pointer) as the output of each component: the Rust passes
//! `&mut perm[(num - 1) as usize..]`.

pub mod level;

pub use level::*;

use crate::HighsInt;

/// HiGHS: rcm.cpp line 9 `const HighsInt offset = 1`
pub const offset: HighsInt = 1;

/// HiGHS: rcm.cpp::Highs_rcm (lines 347-566)
/// Returns 1 on a bad `iccsze` / queue state (the C++ prints and returns 1).
pub fn Highs_rcm(
    root: HighsInt,
    adj_num: HighsInt,
    adj_row: &[HighsInt],
    adj: &[HighsInt],
    mask: &mut [HighsInt],
    perm: &mut [HighsInt],
    iccsze: &mut HighsInt,
    node_num: HighsInt,
) -> HighsInt {
    // HiGHS: rcm.cpp lines 428-441 local declarations
    let mut deg: Vec<HighsInt>;
    let mut fnbr: HighsInt;
    let mut i: HighsInt;
    let mut j: HighsInt;
    let mut jstop: HighsInt;
    let mut jstrt: HighsInt;
    let mut k: HighsInt;
    let mut l: HighsInt;
    let mut lbegin: HighsInt;
    let mut lnbr: HighsInt;
    let mut lperm: HighsInt;
    let mut lvlend: HighsInt;
    let mut nbr: HighsInt;
    let mut node: HighsInt;

    // HiGHS: rcm.cpp lines 443-449
    //  If node_num out of bounds, something is wrong.
    if node_num < 1 {
        // HIGHS_ORDERING_PRINT dropped (lines 446-447)
        return 1;
    }
    // HiGHS: rcm.cpp lines 451-457
    //  If the root is out of bounds, something is wrong.
    if root < 1 || node_num < root {
        // HIGHS_ORDERING_PRINT dropped (lines 454-455)
        return 1;
    }
    // HiGHS: rcm.cpp line 461
    //  Allocate memory for the degree array.
    deg = vec![0; node_num as usize];
    // HiGHS: rcm.cpp lines 463-466
    //  Find the degrees of the nodes in the component specified by MASK and ROOT.
    Highs_rcm_degree(
        root, adj_num, adj_row, adj, mask, &mut deg, iccsze, perm, node_num,
    );
    // HiGHS: rcm.cpp lines 468-475
    //  If the connected component size is less than 1, something is wrong.
    if *iccsze < 1 {
        // HIGHS_ORDERING_PRINT dropped (lines 471-473)
        return 1;
    }
    // HiGHS: rcm.cpp line 479
    //  Set the mask value for the root.
    mask[(root - 1) as usize] = 0;
    // HiGHS: rcm.cpp lines 481-486
    //  If the connected component is a singleton, there is no ordering necessary.
    if *iccsze == 1 {
        // delete[] deg; (line 484) — Drop handles this
        return 0;
    }
    // HiGHS: rcm.cpp lines 488-494
    //  Carry out the reordering.
    //  LBEGIN and LVLEND point to the beginning and
    //  the end of the current level respectively.
    lvlend = 0;
    lnbr = 1;

    // HiGHS: rcm.cpp line 496
    while lvlend < lnbr {
        lbegin = lvlend + 1;
        lvlend = lnbr;

        // HiGHS: rcm.cpp line 500
        i = lbegin;
        while i <= lvlend {
            // HiGHS: rcm.cpp line 504
            //  For each node in the current level...
            node = perm[(i - 1) as usize];
            jstrt = adj_row[(node - 1) as usize] + offset;
            jstop = adj_row[node as usize] - 1 + offset;
            // HiGHS: rcm.cpp lines 508-513
            //  Find the unnumbered neighbors of NODE.
            //  FNBR and LNBR point to the first and last neighbors
            //  of the current node in PERM.
            fnbr = lnbr + 1;

            // HiGHS: rcm.cpp line 515
            j = jstrt;
            while j <= jstop {
                nbr = adj[(j - 1) as usize] + offset;

                if mask[(nbr - 1) as usize] != 0 {
                    lnbr = lnbr + 1;
                    mask[(nbr - 1) as usize] = 0;
                    perm[(lnbr - 1) as usize] = nbr;
                }
                j += 1;
            }
            // HiGHS: rcm.cpp lines 525-529
            //  If no neighbors, skip to next node in this level.
            if lnbr <= fnbr {
                i += 1;
                continue;
            }
            // HiGHS: rcm.cpp lines 531-552
            //  Sort the neighbors of NODE in increasing order by degree.
            //  Linear insertion is used.
            k = fnbr;

            while k < lnbr {
                l = k;
                k = k + 1;
                nbr = perm[(k - 1) as usize];

                while fnbr < l {
                    lperm = perm[(l - 1) as usize];

                    if deg[(lperm - 1) as usize] <= deg[(nbr - 1) as usize] {
                        break;
                    }

                    perm[l as usize] = lperm;
                    l = l - 1;
                }
                perm[l as usize] = nbr;
            }
            i += 1;
        }
    }
    // HiGHS: rcm.cpp lines 556-559
    //  We now have the Cuthill-McKee ordering.
    //  Reverse it to get the Reverse Cuthill-McKee ordering.
    Highs_rcm_i4vec_reverse(*iccsze, perm);
    // HiGHS: rcm.cpp line 563
    //  Free memory.
    // delete[] deg; — Drop handles this

    0
}

/// HiGHS: rcm.cpp::Highs_genrcm (lines 758-879)
/// Output `perm` is 0-based (the last loop subtracts `offset`).
pub fn Highs_genrcm(
    node_num: HighsInt,
    adj_num: HighsInt,
    adj_row: &[HighsInt],
    adj: &[HighsInt],
    perm: &mut [HighsInt],
) -> HighsInt {
    // HiGHS: rcm.cpp lines 821-827 local declarations
    let mut i: HighsInt;
    let mut iccsze: HighsInt = 0;
    let mut level_num: HighsInt = 0;
    let mut level_row: Vec<HighsInt>;
    let mut mask: Vec<HighsInt>;
    let mut num: HighsInt;
    let mut root: HighsInt;

    // HiGHS: rcm.cpp lines 829-830
    level_row = vec![0; (node_num + 1) as usize];
    mask = vec![0; node_num as usize];

    // HiGHS: rcm.cpp lines 832-834
    i = 0;
    while i < node_num {
        mask[i as usize] = 1;
        i += 1;
    }

    // HiGHS: rcm.cpp line 836
    num = 1;

    // HiGHS: rcm.cpp line 838
    i = 0;
    while i < node_num {
        // HiGHS: rcm.cpp lines 839-842
        //  For each masked connected component...
        if mask[i as usize] != 0 {
            root = i + 1;
            // HiGHS: rcm.cpp lines 844-849
            //  Find a pseudo-peripheral node ROOT.  The level structure found by
            //  ROOT_FIND is stored starting at PERM(NUM).
            Highs_rcm_root_find(
                &mut root,
                adj_num,
                adj_row,
                adj,
                &mut mask,
                &mut level_num,
                &mut level_row,
                &mut perm[(num - 1) as usize..],
                node_num,
            );
            // HiGHS: rcm.cpp lines 851-858
            //  RCM orders the component using ROOT as the starting node.
            if Highs_rcm(
                root,
                adj_num,
                adj_row,
                adj,
                &mut mask,
                &mut perm[(num - 1) as usize..],
                &mut iccsze,
                node_num,
            ) != 0
            {
                // delete[] level_row; delete[] mask; (lines 855-856) — Drop handles this
                return 1;
            }

            num = num + iccsze;
            // HiGHS: rcm.cpp lines 862-866
            //  We can stop once every node is in one of the connected components.
            if node_num < num {
                break;
            }
        }
        i += 1;
    }

    // HiGHS: rcm.cpp lines 870-873
    // convert to zero-based numbering
    i = 0;
    while i < node_num {
        perm[i as usize] -= offset;
        i += 1;
    }

    // HiGHS: rcm.cpp lines 875-876
    // delete[] level_row; delete[] mask; — Drop handles this

    0
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn genrcm_path_graph() {
        // path 0-1-2-3 in 0-based CSR (as chooseOrdering passes it)
        let adj_row = [0, 1, 3, 5, 6];
        let adj = [1, 0, 2, 1, 3, 2];
        let mut perm = [0; 4];
        assert_eq!(Highs_genrcm(4, 6, &adj_row, &adj, &mut perm), 0);
        assert!(perm == [0, 1, 2, 3] || perm == [3, 2, 1, 0]);
    }
}
