//! HiGHS source: extern/rcm/rcm.cpp (lines 13-344, 569-755)
//!
//! Translated from RCM as vendored by HiGHS in `extern/rcm` (MIT License, Copyright (c)
//! 2005 Alan George, Joseph Liu, John Burkardt). See THIRD_PARTY_NOTICES.md.
//!
//! Component degrees, rooted level structures and the pseudo-peripheral root
//! search of RCM (all 1-based, see [`super`]).

use super::offset;
use crate::HighsInt;

/// HiGHS: rcm.cpp::Highs_rcm_degree (lines 13-150)
pub fn Highs_rcm_degree(
    root: HighsInt,
    adj_num: HighsInt,
    adj_row: &[HighsInt],
    adj: &[HighsInt],
    mask: &mut [HighsInt],
    deg: &mut [HighsInt],
    iccsze: &mut HighsInt,
    ls: &mut [HighsInt],
    node_num: HighsInt,
) {
    let mut i: HighsInt;
    let mut ideg: HighsInt;
    let mut j: HighsInt;
    let mut jstop: HighsInt;
    let mut jstrt: HighsInt;
    let mut lbegin: HighsInt;
    let mut lvlend: HighsInt;
    let mut lvsize: HighsInt;
    let mut nbr: HighsInt;
    let mut node: HighsInt;

    // The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
    // Use a separate array instead, to keep adj_row const.
    let mut mark = vec![1 as HighsInt; node_num as usize];

    ls[0] = root;
    mark[(root - 1) as usize] = -1;
    lvlend = 0;
    *iccsze = 1;

    // LBEGIN is the pointer to the beginning of the current level, and
    // LVLEND points to the end of this level.
    loop {
        lbegin = lvlend + 1;
        lvlend = *iccsze;

        // Find the degrees of nodes in the current level,
        // and at the same time, generate the next level.
        i = lbegin;
        while i <= lvlend {
            node = ls[(i - 1) as usize];
            jstrt = adj_row[(node - 1) as usize] + offset;
            jstop = adj_row[node as usize] - 1 + offset;
            ideg = 0;

            j = jstrt;
            while j <= jstop {
                nbr = adj[(j - 1) as usize] + offset;

                if mask[(nbr - 1) as usize] != 0 {
                    ideg = ideg + 1;

                    if 0 <= mark[(nbr - 1) as usize] {
                        mark[(nbr - 1) as usize] = -1;
                        *iccsze = *iccsze + 1;
                        ls[(*iccsze - 1) as usize] = nbr;
                    }
                }
                j += 1;
            }
            deg[(node - 1) as usize] = ideg;
            i += 1;
        }

        // Compute the current level width.
        lvsize = *iccsze - lvlend;

        // If the current level width is nonzero, generate another level.
        if lvsize == 0 {
            break;
        }
    }

    // Reset ADJ_ROW to its correct sign and return.
    // (The C++ loop body is a no-op; kept for structural fidelity.)
    i = 0;
    while i < *iccsze {
        let _ = ls[i as usize] - 1;
        i += 1;
    }

    // delete[] mark; (Rust: Drop)
}

/// HiGHS: rcm.cpp::Highs_rcm_level_set (lines 208-344)
pub fn Highs_rcm_level_set(
    root: HighsInt,
    adj_num: HighsInt,
    adj_row: &[HighsInt],
    adj: &[HighsInt],
    mask: &mut [HighsInt],
    level_num: &mut HighsInt,
    level_row: &mut [HighsInt],
    level: &mut [HighsInt],
    node_num: HighsInt,
) {
    let mut i: HighsInt;
    let mut iccsze: HighsInt;
    let mut j: HighsInt;
    let mut jstop: HighsInt;
    let mut jstrt: HighsInt;
    let mut lbegin: HighsInt;
    let mut lvlend: HighsInt;
    let mut lvsize: HighsInt;
    let mut nbr: HighsInt;
    let mut node: HighsInt;

    mask[(root - 1) as usize] = 0;
    level[0] = root;
    *level_num = 0;
    lvlend = 0;
    iccsze = 1;

    // LBEGIN is the pointer to the beginning of the current level, and
    // LVLEND points to the end of this level.
    loop {
        lbegin = lvlend + 1;
        lvlend = iccsze;
        *level_num = *level_num + 1;
        level_row[(*level_num - 1) as usize] = lbegin;

        // Generate the next level by finding all the masked neighbors of nodes
        // in the current level.
        i = lbegin;
        while i <= lvlend {
            node = level[(i - 1) as usize];
            jstrt = adj_row[(node - 1) as usize] + offset;
            jstop = adj_row[node as usize] - 1 + offset;

            j = jstrt;
            while j <= jstop {
                nbr = adj[(j - 1) as usize] + offset;

                if mask[(nbr - 1) as usize] != 0 {
                    iccsze = iccsze + 1;
                    level[(iccsze - 1) as usize] = nbr;
                    mask[(nbr - 1) as usize] = 0;
                }
                j += 1;
            }
            i += 1;
        }

        // Compute the current level width (the number of nodes encountered.)
        // If it is positive, generate the next level.
        lvsize = iccsze - lvlend;

        if lvsize <= 0 {
            break;
        }
    }

    level_row[*level_num as usize] = lvlend + 1;

    // Reset MASK to 1 for the nodes in the level structure.
    i = 0;
    while i < iccsze {
        mask[(level[i as usize] - 1) as usize] = 1;
        i += 1;
    }
}

/// HiGHS: rcm.cpp::Highs_rcm_i4vec_reverse (lines 153-205)
pub fn Highs_rcm_i4vec_reverse(n: HighsInt, a: &mut [HighsInt]) {
    let mut i: HighsInt;
    let mut j: HighsInt;

    i = 0;
    while i < n / 2 {
        j = a[i as usize];
        a[i as usize] = a[(n - 1 - i) as usize];
        a[(n - 1 - i) as usize] = j;
        i += 1;
    }
}

/// HiGHS: rcm.cpp::Highs_rcm_root_find (lines 569-755)
pub fn Highs_rcm_root_find(
    root: &mut HighsInt,
    adj_num: HighsInt,
    adj_row: &[HighsInt],
    adj: &[HighsInt],
    mask: &mut [HighsInt],
    level_num: &mut HighsInt,
    level_row: &mut [HighsInt],
    level: &mut [HighsInt],
    node_num: HighsInt,
) {
    let mut iccsze: HighsInt;
    let mut j: HighsInt;
    let mut jstrt: HighsInt;
    let mut k: HighsInt;
    let mut kstop: HighsInt;
    let mut kstrt: HighsInt;
    let mut level_num2: HighsInt = 0;
    let mut mindeg: HighsInt;
    let mut nabor: HighsInt;
    let mut ndeg: HighsInt;
    let mut node: HighsInt;

    //
    //  Determine the level structure rooted at ROOT.
    //
    Highs_rcm_level_set(
        *root, adj_num, adj_row, adj, mask, level_num, level_row, level, node_num,
    );
    //
    //  Count the number of nodes in this level structure.
    //
    iccsze = level_row[*level_num as usize] - 1;
    //
    //  Extreme case:
    //    A complete graph has a level set of only a single level.
    //    Every node is equally good (or bad).
    //
    if *level_num == 1 {
        return;
    }
    //
    //  Extreme case:
    //    A "line graph" 0--0--0--0--0 has every node in its only level.
    //    By chance, we've stumbled on the ideal root.
    //
    if *level_num == iccsze {
        return;
    }
    //
    //  Pick any node from the last level that has minimum degree
    //  as the starting point to generate a new level set.
    //
    loop {
        mindeg = iccsze;

        jstrt = level_row[(*level_num - 1) as usize];
        *root = level[(jstrt - 1) as usize];

        if jstrt < iccsze {
            j = jstrt;
            while j <= iccsze {
                node = level[(j - 1) as usize];
                ndeg = 0;
                kstrt = adj_row[(node - 1) as usize] + offset;
                kstop = adj_row[node as usize] - 1 + offset;

                k = kstrt;
                while k <= kstop {
                    nabor = adj[(k - 1) as usize] + offset;
                    if 0 < mask[(nabor - 1) as usize] {
                        ndeg = ndeg + 1;
                    }
                    k += 1;
                }

                if ndeg < mindeg {
                    *root = node;
                    mindeg = ndeg;
                }
                j += 1;
            }
        }
        //
        //  Generate the rooted level structure associated with this node.
        //
        Highs_rcm_level_set(
            *root,
            adj_num,
            adj_row,
            adj,
            mask,
            &mut level_num2,
            level_row,
            level,
            node_num,
        );
        //
        //  If the number of levels did not increase, accept the new ROOT.
        //
        if level_num2 <= *level_num {
            break;
        }

        *level_num = level_num2;
        //
        //  In the unlikely case that ROOT is one endpoint of a line graph,
        //  we can exit now.
        //
        if iccsze <= *level_num {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn reverse_in_place() {
        let mut a = [1, 2, 3, 4, 5];
        Highs_rcm_i4vec_reverse(5, &mut a);
        assert_eq!(a, [5, 4, 3, 2, 1]);
    }

    // path 1-2-3-4 (1-based nodes), 0-based CSR as Highs_genrcm receives it
    const PATH_ROW: [HighsInt; 5] = [0, 1, 3, 5, 6];
    const PATH_ADJ: [HighsInt; 6] = [1, 0, 2, 1, 3, 2];

    #[test]
    fn degree_path_all_unmasked() {
        let mut mask = [1, 1, 1, 1];
        let mut deg = [-7; 4];
        let mut iccsze = 0;
        let mut ls = [0; 4];
        Highs_rcm_degree(
            2,
            6,
            &PATH_ROW,
            &PATH_ADJ,
            &mut mask,
            &mut deg,
            &mut iccsze,
            &mut ls,
            4,
        );
        assert_eq!(iccsze, 4);
        assert_eq!(ls, [2, 1, 3, 4]);
        assert_eq!(deg, [1, 2, 2, 1]);
        assert_eq!(mask, [1, 1, 1, 1]);
    }

    #[test]
    fn degree_skips_masked_node() {
        // node 3 has mask 0: neither counted in deg[2-1] nor reached
        let mut mask = [1, 1, 0, 1];
        let mut deg = [-7; 4];
        let mut iccsze = 0;
        let mut ls = [0; 4];
        Highs_rcm_degree(
            2,
            6,
            &PATH_ROW,
            &PATH_ADJ,
            &mut mask,
            &mut deg,
            &mut iccsze,
            &mut ls,
            4,
        );
        assert_eq!(iccsze, 2);
        assert_eq!(ls[..2], [2, 1]);
        assert_eq!(deg, [1, 1, -7, -7]);
        assert_eq!(mask, [1, 1, 0, 1]);
    }

    #[test]
    fn level_set_path_all_unmasked() {
        let mut mask = [1, 1, 1, 1];
        let mut level_num = 0;
        let mut level_row = [-7; 5];
        let mut level = [0; 4];
        Highs_rcm_level_set(
            2,
            6,
            &PATH_ROW,
            &PATH_ADJ,
            &mut mask,
            &mut level_num,
            &mut level_row,
            &mut level,
            4,
        );
        assert_eq!(level_num, 3);
        assert_eq!(level_row, [1, 2, 4, 5, -7]);
        assert_eq!(level, [2, 1, 3, 4]);
        assert_eq!(mask, [1, 1, 1, 1]);
    }

    #[test]
    fn level_set_stops_at_masked_node() {
        // mask is reset to 1 only for nodes in the level structure; node 3 stays 0
        let mut mask = [1, 1, 0, 1];
        let mut level_num = 0;
        let mut level_row = [-7; 5];
        let mut level = [-7; 4];
        Highs_rcm_level_set(
            1,
            6,
            &PATH_ROW,
            &PATH_ADJ,
            &mut mask,
            &mut level_num,
            &mut level_row,
            &mut level,
            4,
        );
        assert_eq!(level_num, 2);
        assert_eq!(level_row, [1, 2, 3, -7, -7]);
        assert_eq!(level, [1, 2, -7, -7]);
        assert_eq!(mask, [1, 1, 0, 1]);
    }
}
