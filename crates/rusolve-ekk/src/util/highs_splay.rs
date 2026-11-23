//! HiGHS source: util/HighsSplay.h
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Top-down splay tree over array-stored nodes (HPresolve keeps one tree per
//! row; node = nonzero position, key = column index).
//!
//! The C++ templates take three lambdas `get_left(node) -> HighsInt&`,
//! `get_right(node) -> HighsInt&`, `get_key(node)`. Every call site in scope
//! (HPresolve link/unlink/findNonzero) uses `ARleft`, `ARright` and `Acol`,
//! so the Rust functions take those arrays directly:
//! `get_left(n)` -> `left[n as usize]`, `get_right(n)` -> `right[n as usize]`,
//! `get_key(n)` -> `keys[n as usize]`.
//!
//! `HighsInt& root` parameters are `&mut HighsInt` (call sites pass
//! `&mut self.rowroot[row as usize]`, a field disjoint from the arrays).
//!
//! Translation notes for the filler:
//! * `highs_splay` keeps two C++ pointers `lright`/`rleft` that point either
//!   at the locals `Nright`/`Nleft` or at a child slot of some node
//!   (`&get_right(root)` / `&get_left(root)`). Model each as
//!   `enum Slot { Local, Left(HighsInt), Right(HighsInt) }` and write through
//!   a small helper (`Local` writes the matching local). Assignment order
//!   must stay exactly as in C++.
//! * The recursive call `highs_splay_unlink(unlinknode, get_right(root), ...)`
//!   passes a reference into `right` while `right` itself is also borrowed.
//!   Copy in / copy out: `let mut sub = right[root as usize];
//!   highs_splay_unlink(unlinknode, &mut sub, left, right, keys);
//!   right[root as usize] = sub;`. This is equivalent because `root` is not
//!   inside its own right subtree, so the recursion never reads or writes
//!   `right[root]` through the array.

use crate::prelude::*;

/// HiGHS: HighsSplay.h::highs_splay (lines 26-76, 51 lines)
///
/// Returns the new root. `key` is compared with `<` / `>` exactly as in C++.
pub fn highs_splay(
    key: HighsInt,
    root: HighsInt,
    left: &mut [HighsInt],
    right: &mut [HighsInt],
    keys: &[HighsInt],
) -> HighsInt {
    if root == -1 {
        return -1;
    }

    let mut nleft: HighsInt = -1;
    let mut nright: HighsInt = -1;
    // `lright = &Nright` is `None`; afterwards `lright = &get_right(root)` is `Some(root)`.
    let mut lright: Option<HighsInt> = None;
    // `rleft = &Nleft` is `None`; afterwards `rleft = &get_left(root)` is `Some(root)`.
    let mut rleft: Option<HighsInt> = None;
    let mut root = root;

    loop {
        if key < keys[root as usize] {
            let left_child = left[root as usize];
            if left_child == -1 {
                break;
            }
            if key < keys[left_child as usize] {
                let y = left_child;
                left[root as usize] = right[y as usize];
                right[y as usize] = root;
                root = y;
                if left[root as usize] == -1 {
                    break;
                }
            }

            // *rleft = root;
            match rleft {
                None => nleft = root,
                Some(n) => left[n as usize] = root,
            }
            // rleft = &get_left(root);
            rleft = Some(root);
            root = left[root as usize];
        } else if key > keys[root as usize] {
            let right_child = right[root as usize];
            if right_child == -1 {
                break;
            }
            if key > keys[right_child as usize] {
                let y = right_child;
                right[root as usize] = left[y as usize];
                left[y as usize] = root;
                root = y;
                if right[root as usize] == -1 {
                    break;
                }
            }

            // *lright = root;
            match lright {
                None => nright = root,
                Some(n) => right[n as usize] = root,
            }
            // lright = &get_right(root);
            lright = Some(root);
            root = right[root as usize];
        } else {
            break;
        }
    }

    // *lright = get_left(root);
    let v = left[root as usize];
    match lright {
        None => nright = v,
        Some(n) => right[n as usize] = v,
    }
    // *rleft = get_right(root);
    let v = right[root as usize];
    match rleft {
        None => nleft = v,
        Some(n) => left[n as usize] = v,
    }
    left[root as usize] = nright;
    right[root as usize] = nleft;

    root
}

/// HiGHS: HighsSplay.h::highs_splay_link (lines 80-104, 25 lines)
///
/// Equal keys go to the right subtree.
pub fn highs_splay_link(
    linknode: HighsInt,
    root: &mut HighsInt,
    left: &mut [HighsInt],
    right: &mut [HighsInt],
    keys: &[HighsInt],
) {
    if *root == -1 {
        left[linknode as usize] = -1;
        right[linknode as usize] = -1;
        *root = linknode;
        return;
    }

    *root = highs_splay(keys[linknode as usize], *root, left, right, keys);

    if keys[linknode as usize] < keys[*root as usize] {
        left[linknode as usize] = left[*root as usize];
        right[linknode as usize] = *root;
        left[*root as usize] = -1;
    } else {
        debug_assert!(keys[linknode as usize] > keys[*root as usize]);
        right[linknode as usize] = right[*root as usize];
        left[linknode as usize] = *root;
        right[*root as usize] = -1;
    }

    *root = linknode;
}

/// HiGHS: HighsSplay.h::highs_splay_unlink (lines 108-133, 26 lines)
pub fn highs_splay_unlink(
    unlinknode: HighsInt,
    root: &mut HighsInt,
    left: &mut [HighsInt],
    right: &mut [HighsInt],
    keys: &[HighsInt],
) {
    debug_assert!(*root != -1);
    *root = highs_splay(keys[unlinknode as usize], *root, left, right, keys);
    debug_assert!(keys[*root as usize] == keys[unlinknode as usize]);

    // Equal keys go to the right subtree, so the node may sit below the root:
    // recurse into right[root] (copy in / copy out, see module docs).
    if *root != unlinknode {
        let mut sub = right[*root as usize];
        highs_splay_unlink(unlinknode, &mut sub, left, right, keys);
        right[*root as usize] = sub;
        return;
    }

    debug_assert!(*root == unlinknode);

    if left[unlinknode as usize] == -1 {
        *root = right[unlinknode as usize];
    } else {
        *root = highs_splay(
            keys[unlinknode as usize],
            left[unlinknode as usize],
            left,
            right,
            keys,
        );
        right[*root as usize] = right[unlinknode as usize];
    }
}
