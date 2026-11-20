//! HiGHS source: ipm/hipo/auxiliary/Auxiliary.h, ipm/hipo/auxiliary/Auxiliary.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Sparse-pattern and elimination-tree helpers shared by FactorHiGHS and the
//! IPM (`Auxiliary.{h,cpp}`), plus the header templates (`counts2Ptr`,
//! `permuteVector`, `permuteVectorInverse`, `freeVector`) and `Clock`.
//! `printTest` (debug file dump) is not ported.

pub mod logger;
pub mod vector_operations;

use super::{Int, Int64};
use std::time::Instant;

/// HiGHS: Auxiliary.cpp::inversePerm (lines 7-15)
pub fn inversePerm(perm: &[Int], iperm: &mut Vec<Int>) {
    // Given the permutation perm, produce the inverse permutation iperm.
    // perm[i] : i-th entry to use in the new order.
    // iperm[i]: where entry i is located in the new order.
    // iperm is pre-sized by the caller (C++ does not resize it).
    for i in 0..perm.len() as Int {
        iperm[perm[i as usize] as usize] = i;
    }
}

/// HiGHS: Auxiliary.cpp::subtreeSize (lines 17-27)
pub fn subtreeSize(parent: &[Int], sizes: &mut Vec<Int>) {
    // Compute sizes of subtrees of the tree given by parent
    let n = parent.len() as Int;
    sizes.clear();
    sizes.resize(n as usize, 1);

    for i in 0..n {
        let k = parent[i as usize];
        if k != -1 {
            sizes[k as usize] += sizes[i as usize];
        }
    }
}

/// HiGHS: Auxiliary.cpp::transpose (lines 29-54)
/// Pattern-only overload.
pub fn transpose(ptr: &[Int], rows: &[Int], ptrT: &mut Vec<Int>, rowsT: &mut Vec<Int>) {
    // Compute the transpose of the matrix and return it in rowsT and ptrT
    let n = ptr.len() as Int - 1;

    // C++: std::vector<Int> work(n);  (value-initialised to 0)
    let mut work = vec![0 as Int; n as usize];

    // count the entries in each row into work
    for i in 0..ptr[ptr.len() - 1] {
        work[rows[i as usize] as usize] += 1;
    }

    // sum row sums to obtain pointers
    // (ptrT and rowsT are pre-sized by the caller; C++ does not resize them)
    counts2Ptr(ptrT.as_mut_slice(), work.as_mut_slice());

    for j in 0..n {
        for el in ptr[j as usize]..ptr[(j + 1) as usize] {
            let i = rows[el as usize];

            // entry (i,j) becomes entry (j,i)
            let pos = work[i as usize];
            work[i as usize] += 1;
            rowsT[pos as usize] = j;
        }
    }
}

/// HiGHS: Auxiliary.cpp::transpose (lines 56-83)
/// Overload with values.
pub fn transpose_val(
    ptr: &[Int],
    rows: &[Int],
    val: &[f64],
    ptrT: &mut Vec<Int>,
    rowsT: &mut Vec<Int>,
    valT: &mut Vec<f64>,
) {
    // Compute the transpose of the matrix and return it in rowsT, ptrT and valT
    let n = ptr.len() as Int - 1;

    // C++: std::vector<Int> work(n);
    let mut work = vec![0 as Int; n as usize];

    // count the entries in each row into work
    for i in 0..ptr[ptr.len() - 1] {
        work[rows[i as usize] as usize] += 1;
    }

    // sum row sums to obtain pointers
    // (ptrT, rowsT and valT are pre-sized by the caller; C++ does not resize them)
    counts2Ptr(ptrT.as_mut_slice(), work.as_mut_slice());

    for j in 0..n {
        for el in ptr[j as usize]..ptr[(j + 1) as usize] {
            let i = rows[el as usize];

            // entry (i,j) becomes entry (j,i)
            let pos = work[i as usize];
            work[i as usize] += 1;
            rowsT[pos as usize] = j;
            valT[pos as usize] = val[el as usize];
        }
    }
}

/// HiGHS: Auxiliary.cpp::permuteSym (lines 85-154)
pub fn permuteSym(
    iperm: &[Int],
    ptr: &mut Vec<Int>,
    rows: &mut Vec<Int>,
    val: &mut Vec<f64>,
    lower: bool,
) {
    // Symmetric permutation of the lower (upper) triangular matrix M based on
    // inverse permutation iperm. The resulting matrix is lower (upper)
    // triangular, regardless of the input matrix.
    let n = ptr.len() as Int - 1;
    let work = &mut vec![0 as Int; n as usize];
    let use_val = !val.is_empty();

    // go through the columns to count the nonzeros
    for j in 0..n {
        // get new index of column
        let col = iperm[j as usize];

        // go through elements of column
        for el in ptr[j as usize]..ptr[(j + 1) as usize] {
            let i = rows[el as usize];

            // ignore potential entries in upper(lower) triangular part
            if (lower && i < j) || (!lower && i > j) {
                continue;
            }

            // get new index of row
            let row = iperm[i as usize];

            // if only lower triangular part is used, col is smaller than row
            // C++: std::min / std::max -> cpp_min / cpp_max
            let actual_col = if lower {
                crate::util::highs_utils::cpp_min(row, col)
            } else {
                crate::util::highs_utils::cpp_max(row, col)
            };
            work[actual_col as usize] += 1;
        }
    }

    let mut new_ptr = vec![0 as Int; (n + 1) as usize];

    // get column pointers by summing the count of nonzeros in each column.
    // copy column pointers into work
    counts2Ptr(new_ptr.as_mut_slice(), work.as_mut_slice());

    let nnz = new_ptr[new_ptr.len() - 1] as usize;
    let mut new_rows = vec![0 as Int; nnz];
    let mut new_val: Vec<f64> = Vec::new();
    if use_val {
        new_val.resize(nnz, 0.0);
    }

    // go through the columns to assign row indices
    for j in 0..n {
        // get new index of column
        let col = iperm[j as usize];

        // go through elements of column
        for el in ptr[j as usize]..ptr[(j + 1) as usize] {
            let i = rows[el as usize];

            // ignore potential entries in upper triangular part
            if (lower && i < j) || (!lower && i > j) {
                continue;
            }

            // get new index of row
            let row = iperm[i as usize];

            // if only lower triangular part is used, col is smaller than row
            let actual_col = if lower {
                crate::util::highs_utils::cpp_min(row, col)
            } else {
                crate::util::highs_utils::cpp_max(row, col)
            };
            let actual_row = if lower {
                crate::util::highs_utils::cpp_max(row, col)
            } else {
                crate::util::highs_utils::cpp_min(row, col)
            };

            let pos = work[actual_col as usize];
            work[actual_col as usize] += 1;
            new_rows[pos as usize] = actual_row;
            if use_val {
                new_val[pos as usize] = val[el as usize];
            }
        }
    }

    *ptr = new_ptr;
    *rows = new_rows;
    if use_val {
        *val = new_val;
    }
}

/// HiGHS: Auxiliary.cpp::childrenLinkedList (lines 156-173)
pub fn childrenLinkedList(parent: &[Int], head: &mut Vec<Int>, next: &mut Vec<Int>) {
    // Create linked lists of children in elimination tree.
    // parent gives the dependencies of the tree,
    // head[node] is the first child of node,
    // next[head[node]] is the second child,
    // next[next[head[node]]] is the third child...
    // until -1 is reached.
    let n = parent.len() as Int;
    head.clear();
    head.resize(n as usize, -1);
    next.clear();
    next.resize(n as usize, -1);
    let mut node = n - 1;
    while node >= 0 {
        if parent[node as usize] == -1 {
            node -= 1;
            continue;
        }
        next[node as usize] = head[parent[node as usize] as usize];
        head[parent[node as usize] as usize] = node;
        node -= 1;
    }
}

/// HiGHS: Auxiliary.cpp::reverseLinkedList (lines 175-196)
pub fn reverseLinkedList(head: &mut Vec<Int>, next: &mut Vec<Int>) {
    // Reverse the linked list of children of each node.
    // If a node has children (a -> b -> c -> -1), the reverse list contains
    // children (c -> b -> a -> -1).
    let n = head.len() as Int;

    for node in 0..n {
        let mut prev_node: Int = -1;
        let mut curr_node = head[node as usize];
        let mut next_node: Int;

        while curr_node != -1 {
            next_node = next[curr_node as usize];
            next[curr_node as usize] = prev_node;
            prev_node = curr_node;
            curr_node = next_node;
        }

        head[node as usize] = prev_node;
    }
}

/// HiGHS: Auxiliary.cpp::dfsPostorder (lines 198-223)
pub fn dfsPostorder(
    node: Int,
    start: &mut Int,
    head: &mut Vec<Int>,
    next: &[Int],
    order: &mut Vec<Int>,
) {
    // Perform depth first search starting from root node and order the nodes
    // starting from the value start. head and next contain the linked list of
    // children.
    let mut stack: Vec<Int> = Vec::new();
    stack.push(node);

    while !stack.is_empty() {
        let current = *stack.last().unwrap();
        let child = head[current as usize];

        if child == -1 {
            // no children left to order,
            // remove from the stack and order
            stack.pop();
            order[*start as usize] = current;
            *start += 1;
        } else {
            // at least one child left to order,
            // add it to the stack and remove it from the list of children
            stack.push(child);
            head[current as usize] = next[child as usize];
        }
    }
}

/// HiGHS: Auxiliary.cpp::processEdge (lines 225-265)
pub fn processEdge(
    j: Int,
    i: Int,
    first: &[Int],
    maxfirst: &mut Vec<Int>,
    delta: &mut Vec<Int>,
    prevleaf: &mut Vec<Int>,
    ancestor: &mut Vec<Int>,
) {
    // Process edge of skeleton matrix.
    // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

    // j not a leaf of ith row subtree
    if i <= j || first[j as usize] <= maxfirst[i as usize] {
        return;
    }

    // max first[j] so far
    maxfirst[i as usize] = first[j as usize];

    // previous leaf of ith row subtree
    let jprev = prevleaf[i as usize];

    // A(i,j) is in the skeleton matrix
    delta[j as usize] += 1;

    if jprev != -1 {
        // find least common ancestor of jprev and j
        let mut q = jprev;
        while q != ancestor[q as usize] {
            q = ancestor[q as usize];
        }

        // path compression
        let mut sparent: Int;
        let mut s = jprev;
        while s != q {
            sparent = ancestor[s as usize];
            ancestor[s as usize] = q;
            s = sparent;
        }

        // consider overlap
        delta[q as usize] -= 1;
    }

    // previous leaf of ith subtree set to j
    prevleaf[i as usize] = j;
}

/// HiGHS: Auxiliary.cpp::getDiagStart (lines 267-280)
/// HiGHS 1.16 (#3448): Auxiliary.cpp::getDiagStart (lines 275-290; lines 278-279)
/// C++ default arguments: `triang = false`.
/// The `n_blocks == 0` guard is not behind an option: 1.15.1 reads
/// `start.back()` on an empty vector there (UB), so there is no 1.15.1 result to keep.
pub fn getDiagStart(
    n: Int,
    k: Int,
    nb: Int,
    n_blocks: Int,
    start: &mut Vec<Int64>,
    triang: bool,
) -> Int64 {
    // start position of diagonal blocks for blocked dense formats
    if n_blocks == 0 {
        return 0;
    }
    start.clear();
    start.resize(n_blocks as usize, 0);
    for i in 1..n_blocks {
        start[i as usize] = start[(i - 1) as usize] + (nb * (n - (i - 1) * nb)) as Int64;
        if triang {
            start[i as usize] -= (nb * (nb - 1) / 2) as Int64;
        }
    }

    let jb = crate::util::highs_utils::cpp_min(nb, k - (n_blocks - 1) * nb);
    let mut result = start[start.len() - 1] + (n - (n_blocks - 1) * nb) as Int64 * (jb as Int64);
    if triang {
        result -= (jb * (jb - 1) / 2) as Int64;
    }
    result
}

/// HiGHS: Auxiliary.cpp::maxDepthTree (lines 282-301)
pub fn maxDepthTree(parent: &[Int]) -> Int {
    let mut max_depth: Int = 0;
    let n = parent.len() as Int;
    let mut depth = vec![-1 as Int; n as usize];
    for i in 0..n {
        let mut node = i;
        let mut value: Int = 1;
        while node != -1 {
            if value > depth[node as usize] {
                depth[node as usize] = value;
            } else {
                break;
            }

            value += 1;
            node = parent[node as usize];
        }
        if parent[i as usize] == -1 {
            max_depth = crate::util::highs_utils::cpp_max(max_depth, depth[i as usize]);
        }
    }
    max_depth
}

/// HiGHS: Auxiliary.cpp::fullFromLower (lines 303-334)
pub fn fullFromLower(ptrL: &[Int], rowsL: &[Int], ptrF: &mut Vec<Int>, rowsF: &mut Vec<Int>) {
    // Given a sparse matrix in lower triangular format, build the same matrix in
    // full format, without diagonal entries.
    let mut rowsU = vec![0 as Int; rowsL.len()];
    let mut ptrU = vec![0 as Int; ptrL.len()];
    transpose(ptrL, rowsL, &mut ptrU, &mut rowsU);

    let n = ptrL.len() as Int - 1;
    let mut work = vec![0 as Int; n as usize];
    for j in 0..n {
        for el in ptrU[j as usize]..ptrU[(j + 1) as usize] {
            let i = rowsU[el as usize];
            if i == j {
                continue;
            }
            work[j as usize] += 1;
            work[i as usize] += 1;
        }
    }

    ptrF.clear();
    ptrF.resize((n + 1) as usize, 0);
    counts2Ptr(ptrF.as_mut_slice(), work.as_mut_slice());
    rowsF.clear();
    rowsF.resize(ptrF[ptrF.len() - 1] as usize, 0);
    for j in 0..n {
        for el in ptrU[j as usize]..ptrU[(j + 1) as usize] {
            let i = rowsU[el as usize];
            if i == j {
                continue;
            }
            let pos_j = work[j as usize];
            work[j as usize] += 1;
            rowsF[pos_j as usize] = i;
            let pos_i = work[i as usize];
            work[i as usize] += 1;
            rowsF[pos_i as usize] = j;
        }
    }
}

/// HiGHS: Auxiliary.cpp::snFlops (lines 336-340)
pub fn snFlops(size: f64, clique_size: f64) -> f64 {
    (size + clique_size) * (size + clique_size) * size - (size + clique_size) * size * (size + 1.0)
        + size * (size + 1.0) * (2.0 * size + 1.0) / 6.0
}

/// HiGHS: Auxiliary.cpp::snSpops (lines 341-343)
pub fn snSpops(clique_size: f64) -> f64 {
    clique_size * (clique_size + 1.0) / 2.0
}

/// HiGHS: Auxiliary.h::counts2Ptr (template, lines 47-56)
pub fn counts2Ptr<T>(ptr: &mut [T], w: &mut [T])
where
    T: Copy + Default + std::ops::AddAssign + TryFrom<usize> + TryInto<usize>,
{
    let mut temp_nz = T::default();
    let n = w.len();
    for j in 0..n {
        ptr[j] = temp_nz;
        temp_nz += w[j];
        w[j] = ptr[j];
    }
    ptr[n] = temp_nz;
}

/// HiGHS: Auxiliary.h::permuteVector (template, lines 58-62)
pub fn permuteVector<T: Clone>(v: &mut [T], perm: &[Int]) {
    let temp_v = v.to_vec();
    for i in 0..v.len() {
        v[i] = temp_v[perm[i] as usize].clone();
    }
}

/// HiGHS: Auxiliary.h::permuteVectorInverse (template, lines 64-68)
pub fn permuteVectorInverse<T: Clone>(v: &mut [T], iperm: &[Int]) {
    let temp_v = v.to_vec();
    for i in 0..v.len() {
        v[iperm[i] as usize] = temp_v[i].clone();
    }
}

/// HiGHS: Auxiliary.h::freeVector (template, lines 82-86). `clear()` +
/// `shrink_to_fit()`: the capacity matters only for memory reporting.
pub fn freeVector<T>(v: &mut Vec<T>) {
    v.clear();
    v.shrink_to_fit();
}

/// HiGHS: Auxiliary.h `class Clock` (wall clock; only feeds timing output).
#[derive(Clone, Copy, Debug)]
pub struct Clock {
    pub t0: Instant,
}

impl Default for Clock {
    fn default() -> Self {
        Self::new()
    }
}

impl Clock {
    /// HiGHS: Auxiliary.cpp::Clock::Clock (line 345)
    pub fn new() -> Self {
        Self { t0: Instant::now() }
    }

    /// HiGHS: Auxiliary.cpp::Clock::start (line 346)
    pub fn start(&mut self) {
        self.t0 = Instant::now();
    }

    /// HiGHS: Auxiliary.cpp::Clock::stop (lines 347-351)
    pub fn stop(&self) -> f64 {
        self.t0.elapsed().as_secs_f64()
    }
}

#[cfg(test)]
mod tests {
    //! Unit tests.
    use super::*;

    #[test]
    fn inverse_perm_and_transpose() {
        // The C++ callers pre-size iperm / ptrT / rowsT (e.g. Analyse.cpp lines 36-37,
        // 59), and inversePerm / transpose do not resize them.
        let mut iperm = vec![0; 3];
        inversePerm(&[2, 0, 1], &mut iperm);
        assert_eq!(iperm, vec![1, 2, 0]);
        // 2x2 lower pattern: col 0 rows {0,1}, col 1 rows {1}
        let (mut pt, mut rt) = (vec![0; 3], vec![0; 3]);
        transpose(&[0, 2, 3], &[0, 1, 1], &mut pt, &mut rt);
        assert_eq!(pt, vec![0, 1, 3]);
        assert_eq!(rt, vec![0, 0, 1]);
    }

    #[test]
    fn tree_helpers() {
        // chain 0 -> 1 -> 2 (root)
        assert_eq!(maxDepthTree(&[1, 2, -1]), 3);
        let mut sizes = Vec::new();
        subtreeSize(&[1, 2, -1], &mut sizes);
        assert_eq!(sizes, vec![1, 2, 3]);
        assert_eq!(snSpops(2.0), 3.0);
    }

    #[test]
    fn permute_full_and_dfs() {
        // lower 2x2: col 0 rows {0,1} vals {1,2}, col 1 row {1} val {3}; swap 0 <-> 1
        let (mut p, mut r, mut v) = (vec![0, 2, 3], vec![0, 1, 1], vec![1.0, 2.0, 3.0]);
        permuteSym(&[1, 0], &mut p, &mut r, &mut v, true);
        assert_eq!(
            (p, r, v),
            (vec![0, 2, 3], vec![1, 0, 1], vec![2.0, 3.0, 1.0])
        );
        let (mut pf, mut rf) = (Vec::new(), Vec::new());
        fullFromLower(&[0, 2, 3], &[0, 1, 1], &mut pf, &mut rf);
        assert_eq!((pf, rf), (vec![0, 1, 2], vec![1, 0]));
        // parent: 0 -> 2, 1 -> 2
        let (mut head, mut next) = (Vec::new(), Vec::new());
        childrenLinkedList(&[2, 2, -1], &mut head, &mut next);
        assert_eq!((&head, &next), (&vec![-1, -1, 0], &vec![1, -1, -1]));
        let (mut h2, mut n2) = (head.clone(), next.clone());
        reverseLinkedList(&mut h2, &mut n2);
        assert_eq!((h2, n2), (vec![-1, -1, 1], vec![-1, 0, -1]));
        let (mut start, mut order) = (0, vec![0; 3]);
        dfsPostorder(2, &mut start, &mut head, &next, &mut order);
        assert_eq!((start, order), (3, vec![0, 1, 2]));
        // 5x5 dense, blocks of 2: full 17 entries, lower triangle 15
        let mut st = Vec::new();
        assert_eq!(getDiagStart(5, 5, 2, 3, &mut st, false), 17);
        assert_eq!(st, vec![0, 10, 16]);
        assert_eq!(getDiagStart(5, 5, 2, 3, &mut st, true), 15);
        assert_eq!(st, vec![0, 9, 14]);
    }

    #[test]
    fn h16_hipo_diag_start_zero_blocks() {
        // HiGHS 1.16 #3448: Auxiliary.cpp lines 278-279 return 0 before start.assign when n_blocks == 0.
        // The start vector is not touched on that path.
        let mut st = vec![7, 8];
        assert_eq!(getDiagStart(5, 5, 2, 0, &mut st, false), 0);
        assert_eq!(getDiagStart(5, 5, 2, 0, &mut st, true), 0);
        assert_eq!(st, vec![7, 8]);
    }
}
