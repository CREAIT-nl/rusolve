//! HiGHS source: ipm/hipo/factorhighs/Analyse.cpp (lines 948-1302)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Children reordering, clique blocks, overflow check, stack size and the
//! `run` driver that fills the `Symbolic`.

use super::super::super::{Int, Int64};
use super::super::symbolic::Symbolic;
use super::super::*;
use super::Analyse;
use crate::ipm::hipo::auxiliary::{
    childrenLinkedList, dfsPostorder, getDiagStart, inversePerm, maxDepthTree, permuteVector,
    transpose,
};
use crate::util::highs_utils::{cpp_max, cpp_min};

impl<'a> Analyse<'a> {
    /// HiGHS: Analyse.cpp::Analyse::reorderChildren (lines 948-1108)
    pub fn reorderChildren(&mut self) {
        let sn_count = self.sn_count_;
        let mut clique_entries = vec![0i64; sn_count as usize];
        let mut frontal_entries = vec![0i64; sn_count as usize];
        let mut storage = vec![0f64; sn_count as usize];
        let mut storage_factors = vec![0f64; sn_count as usize];

        // initialise data of supernodes
        for sn in 0..sn_count {
            // supernode size
            let sz = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];

            // frontal size
            let fr = self.col_count_[self.sn_start_[sn as usize] as usize];

            // compute storage based on format used
            let mut fe = 0i64;
            let mut ce = 0i64;
            self.computeStorage(fr, sz, &mut fe, &mut ce);
            frontal_entries[sn as usize] = fe;
            clique_entries[sn as usize] = ce;

            // compute number of entries in factors within the subtree
            storage_factors[sn as usize] += frontal_entries[sn as usize] as f64;
            if self.sn_parent_[sn as usize] != -1 {
                let parent = self.sn_parent_[sn as usize] as usize;
                storage_factors[parent] += storage_factors[sn as usize];
            }
        }

        // linked lists of children
        let mut head: Vec<Int> = Vec::new();
        let mut next: Vec<Int> = Vec::new();
        childrenLinkedList(&self.sn_parent_, &mut head, &mut next);

        // go through the supernodes
        for sn in 0..sn_count {
            // leaf node
            if head[sn as usize] == -1 {
                storage[sn as usize] =
                    (frontal_entries[sn as usize] + clique_entries[sn as usize]) as f64;
                continue;
            }

            // save children and values to sort
            let mut children: Vec<(Int, f64)> = Vec::new();
            let mut child = head[sn as usize];
            while child != -1 {
                let value = storage[child as usize]
                    - clique_entries[child as usize] as f64
                    - storage_factors[child as usize];
                children.push((child, value));
                child = next[child as usize];
            }

            // sort children in decreasing order of the values
            // PORT NOTE: std::sort, keys can tie (siblings with equal values).
            // The tie order is part of the result (it fixes the postorder), so
            // the oracle's std::sort (libstdc++ of GCC 13, the fixture build) is
            // reproduced by `libstdcxx_sort` below. A stable sort_by matches it
            // only up to 16 children.
            libstdcxx_sort(&mut children, |a: &(Int, f64), b: &(Int, f64)| a.1 > b.1);

            // modify linked lists with new order of children
            head[sn as usize] = children[0].0;
            for i in 0..(children.len() as Int - 1) {
                next[children[i as usize].0 as usize] = children[(i + 1) as usize].0;
            }
            next[children[children.len() - 1].0 as usize] = -1;
        }

        // =================================================
        // Create supernodal permutation
        // =================================================
        // build supernodal permutation with dfs
        let mut sn_perm = vec![0i32; sn_count as usize];
        let mut start: Int = 0;
        for sn in 0..sn_count {
            if self.sn_parent_[sn as usize] == -1 {
                dfsPostorder(sn, &mut start, &mut head, &next, &mut sn_perm);
            }
        }

        // =================================================
        // Create nodal permutation
        // =================================================
        // Given the supernodal permutation, find the nodal permutation

        // permutation to apply to the existing one
        let mut new_perm = vec![0i32; self.n_ as usize];

        // index to write into new_perm
        start = 0;

        for i in 0..sn_count {
            let sn = sn_perm[i as usize];
            for j in self.sn_start_[sn as usize]..self.sn_start_[(sn + 1) as usize] {
                new_perm[start as usize] = j;
                start += 1;
            }
        }

        // obtain inverse permutation
        let mut new_iperm = vec![0i32; self.n_ as usize];
        inversePerm(&new_perm, &mut new_iperm);

        // =================================================
        // Create new sn elimination tree
        // =================================================
        let mut isn_perm = vec![0i32; sn_count as usize];
        inversePerm(&sn_perm, &mut isn_perm);
        let mut new_sn_parent = vec![0i32; sn_count as usize];
        for i in 0..sn_count {
            if self.sn_parent_[i as usize] != -1 {
                new_sn_parent[isn_perm[i as usize] as usize] =
                    isn_perm[self.sn_parent_[i as usize] as usize];
            } else {
                new_sn_parent[isn_perm[i as usize] as usize] = -1;
            }
        }

        // =================================================
        // Create new snBelong
        // =================================================

        // build new snBelong, i.e., the sn to which each colum belongs
        for sn in 0..sn_count {
            for i in self.sn_start_[sn as usize]..self.sn_start_[(sn + 1) as usize] {
                self.sn_belong_[i as usize] = isn_perm[sn as usize];
            }
        }
        permuteVector(&mut self.sn_belong_, &new_perm);

        // permute other vectors that may be needed
        permuteVector(&mut self.col_count_, &new_perm);
        permuteVector(&mut self.sn_indices_, &sn_perm);

        // =================================================
        // Create new snStart
        // =================================================
        let mut cols_per_sn = vec![0i32; sn_count as usize];

        // compute size of each supernode
        for sn in 0..sn_count {
            cols_per_sn[sn as usize] =
                self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];
        }

        // permute according to new order of supernodes
        permuteVector(&mut cols_per_sn, &sn_perm);

        // sum number of columns to obtain pointers
        for i in 0..(sn_count - 1) {
            cols_per_sn[(i + 1) as usize] += cols_per_sn[i as usize];
        }

        for i in 0..sn_count {
            self.sn_start_[(i + 1) as usize] = cols_per_sn[i as usize];
        }

        // =================================================
        // Save new data
        // =================================================

        // Overwrite previous data
        self.sn_parent_ = new_sn_parent;

        // Permute matrix based on new permutation
        self.permute(&new_iperm);

        // double transpose to sort columns and compute lower part
        transpose(
            &self.ptr_upper_,
            &self.rows_upper_,
            &mut self.ptr_lower_,
            &mut self.rows_lower_,
        );
        transpose(
            &self.ptr_lower_,
            &self.rows_lower_,
            &mut self.ptr_upper_,
            &mut self.rows_upper_,
        );

        // Update perm and iperm
        permuteVector(&mut self.perm_, &new_perm);
        inversePerm(&self.perm_, &mut self.iperm_);
    }

    /// HiGHS: Analyse.cpp::Analyse::computeBlockStart (lines 1110-1124)
    pub fn computeBlockStart(&mut self) {
        self.clique_block_start_
            .resize(self.sn_count_ as usize, Vec::new());
        // compute starting position of each block of columns in the clique, for
        // each supernode
        for sn in 0..self.sn_count_ {
            let sn_size = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];
            let ldf = (self.ptr_sn_[(sn + 1) as usize] - self.ptr_sn_[sn as usize]) as Int;
            let ldc = ldf - sn_size;
            let n_blocks = (ldc - 1) / self.nb_ + 1;

            let mut start_vec: Vec<Int64> = Vec::new();
            let schur_size = getDiagStart(ldc, ldc, self.nb_, n_blocks, &mut start_vec, false);
            start_vec.push(schur_size);
            self.clique_block_start_[sn as usize] = start_vec;
        }
    }

    /// HiGHS: Analyse.cpp::Analyse::checkOverflow (lines 1126-1143)
    pub fn checkOverflow(&self) -> Int {
        // In order to use 32-bit BLAS, any data accessed by BLAS must be addressable
        // using 32-bit integer offset. If BLAS is given a pointer double* A, the
        // distance between the first and last entry of A used by BLAS needs to be
        // smaller than int32_limit. Since the matrices are stored in blocked data
        // structures, and BLAS only uses contiguous data from a given block of
        // columns, we need to impose that:
        //   front_size * min(block_size, sn_size) <= int32_limit

        // Analyse.cpp:17 `const Int64 int32_limit = std::numeric_limits<int32_t>::max();`
        const int32_limit: Int64 = i32::MAX as Int64;

        for sn in 0..self.sn_count_ {
            let sn_size = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];
            let front_size = (self.ptr_sn_[(sn + 1) as usize] - self.ptr_sn_[sn as usize]) as Int;

            if (front_size as Int64) * cpp_min(sn_size, self.nb_) as Int64 > int32_limit {
                return 1;
            }
        }

        0
    }

    /// HiGHS: Analyse.cpp::Analyse::computeStackSize (lines 1145-1209)
    pub fn computeStackSize(&mut self) {
        // Compute the minimum size of the stack to process the elimination tree
        // serially.

        let sn_count = self.sn_count_;
        let mut clique_entries = vec![0i64; sn_count as usize];
        let mut stack_subtrees = vec![0i64; sn_count as usize];
        let mut total_frontal: Int64 = 0;

        // initialise data of supernodes
        for sn in 0..sn_count {
            // supernode size
            let sz = self.sn_start_[(sn + 1) as usize] - self.sn_start_[sn as usize];

            // frontal size
            let fr = (self.ptr_sn_[(sn + 1) as usize] - self.ptr_sn_[sn as usize]) as Int;

            let mut frontal_entries: Int64 = 0;
            let mut ce: Int64 = 0;

            // compute storage based on format used
            self.computeStorage(fr, sz, &mut frontal_entries, &mut ce);
            clique_entries[sn as usize] = ce;

            total_frontal += frontal_entries;
        }

        // linked lists of children
        let mut head: Vec<Int> = Vec::new();
        let mut next: Vec<Int> = Vec::new();
        childrenLinkedList(&self.sn_parent_, &mut head, &mut next);

        // go through the supernodes
        for sn in 0..sn_count {
            // leaf node
            if head[sn as usize] == -1 {
                stack_subtrees[sn as usize] = clique_entries[sn as usize];
                continue;
            }

            // Compute storage
            // storage is found as max(storage_1,storage_2), where
            // storage_1 = max_j stack_size[j] + \sum_{k up to j-1} clique_entries[k]
            // storage_2 = clique_total_entries (including node itself)

            let mut clique_partial_entries: Int64 = 0;
            let mut storage_1: Int64 = 0;

            let mut child = head[sn as usize];
            while child != -1 {
                let current = stack_subtrees[child as usize] + clique_partial_entries;

                clique_partial_entries += clique_entries[child as usize];
                storage_1 = cpp_max(storage_1, current);

                child = next[child as usize];
            }

            let storage_2 = clique_partial_entries + clique_entries[sn as usize];

            stack_subtrees[sn as usize] = cpp_max(storage_1, storage_2);
            self.max_stack_size_ = cpp_max(self.max_stack_size_, stack_subtrees[sn as usize]);
        }

        // minimum storage in serial is equal to the space needed to store the
        // factorisation and the maximum size of the stack. Times 8 to obtain the
        // number of bytes.
        self.serial_storage_ = ((total_frontal + self.max_stack_size_) * 8) as f64;
    }

    /// HiGHS: Analyse.cpp::Analyse::run (lines 1211-1302)
    pub fn run(&mut self, S: &mut Symbolic) -> Int {
        // Perform analyse phase and store the result into the symbolic object S.
        // After Run returns, the Analyse object is not valid.

        if !self.ready_ {
            return kRetGeneric;
        }

        // HIPO_CLOCK_CREATE; (line 1217) — level-0 no-op, dropped.
        // HIPO_CLOCK_START(2); (line 1219) — level-0 no-op, dropped.
        self.permute(&self.iperm_.clone());
        self.eTree();
        self.postorder();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalyseTree); (line 1223) — level-0 no-op, dropped.

        // HIPO_CLOCK_START(2); (line 1225) — level-0 no-op, dropped.
        self.colCount();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalyseCount); (line 1227) — level-0 no-op, dropped.

        // HIPO_CLOCK_START(2); (line 1229) — level-0 no-op, dropped.
        self.fundamentalSupernodes();
        self.relaxSupernodes();
        self.afterRelaxSn();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalyseSn); (line 1233) — level-0 no-op, dropped.

        // HIPO_CLOCK_START(2); (line 1235) — level-0 no-op, dropped.
        self.reorderChildren();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalyseReorder); (line 1237) — level-0 no-op, dropped.

        // HIPO_CLOCK_START(2); (line 1239) — level-0 no-op, dropped.
        self.snPattern();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalysePattern); (line 1241) — level-0 no-op, dropped.

        // HIPO_CLOCK_START(2); (line 1243) — level-0 no-op, dropped.
        self.relativeIndCols();
        self.relativeIndClique();
        // HIPO_CLOCK_STOP(2, data_, kTimeAnalyseRelInd); (line 1246) — level-0 no-op, dropped.

        self.computeBlockStart();
        self.computeCriticalPath();
        self.computeStackSize();

        // move relevant stuff into S
        S.n_ = self.n_;
        S.sn_ = self.sn_count_;
        S.nz_ = self.nz_factor_;
        S.fillin_ = self.nz_factor_ as f64 / self.nz_ as f64;
        S.artificial_nz_ = self.artificial_nz_;
        S.artificial_ops_ = self.dense_ops_ - self.dense_ops_norelax_;
        S.spops_ = self.sparse_ops_;
        S.critops_ = self.critical_ops_;
        // std::max_element returns first extreme; write the loop
        let mut largest_front: Int = self.sn_indices_[0];
        for i in 1..self.sn_indices_.len() {
            if self.sn_indices_[i] > largest_front {
                largest_front = self.sn_indices_[i];
            }
        }
        S.largest_front_ = largest_front;
        S.serial_storage_ = self.serial_storage_;
        S.flops_ = self.dense_ops_;
        S.block_size_ = self.nb_;
        S.max_stack_size_ = self.max_stack_size_;
        S.tree_depth_ = maxDepthTree(&self.sn_parent_);

        // compute largest supernode
        let mut sn_size: Vec<Int> = self.sn_start_[1..].to_vec();
        for i in (1..self.sn_count_).rev() {
            sn_size[i as usize] -= sn_size[(i - 1) as usize];
        }
        // std::max_element returns first extreme; write the loop
        let mut largest_sn: Int = sn_size[0];
        for i in 1..sn_size.len() {
            if sn_size[i] > largest_sn {
                largest_sn = sn_size[i];
            }
        }
        S.largest_sn_ = largest_sn;

        // build statistics about supernodes size
        for &i in &sn_size {
            if i == 1 {
                S.sn_size_1_ += 1;
            }
            if i <= 10 {
                S.sn_size_10_ += 1;
            }
            if i <= 100 {
                S.sn_size_100_ += 1;
            }
        }

        if self.checkOverflow() != 0 {
            if let Some(logger) = self.logger_ {
                logger.printe("Integer overflow in analyse phase\n");
            }
            return kRetIntOverflow;
        }

        // permute signs of pivots
        S.pivot_sign_ = std::mem::take(&mut self.signs_);
        permuteVector(&mut S.pivot_sign_, &self.perm_);

        S.iperm_ = std::mem::take(&mut self.iperm_);
        S.rows_ = std::mem::take(&mut self.rows_sn_);
        S.ptr_ = std::mem::take(&mut self.ptr_sn_);
        S.sn_parent_ = std::mem::take(&mut self.sn_parent_);
        S.sn_start_ = std::mem::take(&mut self.sn_start_);
        S.relind_cols_ = std::mem::take(&mut self.relind_cols_);
        S.relind_clique_ = std::mem::take(&mut self.relind_clique_);
        S.consecutive_sums_ = std::mem::take(&mut self.consecutive_sums_);
        S.clique_block_start_ = std::mem::take(&mut self.clique_block_start_);

        // HIPO_CLOCK_STOP(1, data_, kTimeAnalyse); (line 1299) — level-0 no-op, dropped.

        kRetOk
    }
}

// =====================================================================
// libstdc++ std::sort (GCC 13, bits/stl_algo.h + bits/stl_heap.h), used by
// reorderChildren. `comp(a, b)` is the C++ comparator ("a before b"). The
// element type is Copy, so each _GLIBCXX_MOVE is a copy. Indices are into the
// whole slice; [first, last) as in the C++ iterators.
// =====================================================================

/// std::__sort
fn libstdcxx_sort<T: Copy, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut comp: F) {
    let (first, last) = (0usize, v.len());
    if first != last {
        // std::__lg(n) * 2, __lg(n) = floor(log2(n))
        let n = last - first;
        let lg = (usize::BITS - 1 - n.leading_zeros()) as usize;
        stdcxx_introsort_loop(v, first, last, lg * 2, &mut comp);
        stdcxx_final_insertion_sort(v, first, last, &mut comp);
    }
}

const STDCXX_S_THRESHOLD: usize = 16;

/// std::__introsort_loop
fn stdcxx_introsort_loop<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    comp: &mut F,
) {
    while last - first > STDCXX_S_THRESHOLD {
        if depth_limit == 0 {
            // std::__partial_sort(first, last, last): __heap_select + __sort_heap
            stdcxx_make_heap(v, first, last, comp);
            stdcxx_sort_heap(v, first, last, comp);
            return;
        }
        depth_limit -= 1;
        let cut = stdcxx_unguarded_partition_pivot(v, first, last, comp);
        stdcxx_introsort_loop(v, cut, last, depth_limit, comp);
        last = cut;
    }
}

/// std::__move_median_to_first
fn stdcxx_move_median_to_first<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    comp: &mut F,
) {
    if comp(&v[a], &v[b]) {
        if comp(&v[b], &v[c]) {
            v.swap(result, b);
        } else if comp(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if comp(&v[a], &v[c]) {
        v.swap(result, a);
    } else if comp(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

/// std::__unguarded_partition
fn stdcxx_unguarded_partition<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    comp: &mut F,
) -> usize {
    loop {
        while comp(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while comp(&v[pivot], &v[last]) {
            last -= 1;
        }
        if !(first < last) {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

/// std::__unguarded_partition_pivot
fn stdcxx_unguarded_partition_pivot<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) -> usize {
    let mid = first + (last - first) / 2;
    stdcxx_move_median_to_first(v, first, first + 1, mid, last - 1, comp);
    stdcxx_unguarded_partition(v, first + 1, last, first, comp)
}

/// std::__unguarded_linear_insert
fn stdcxx_unguarded_linear_insert<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    mut last: usize,
    comp: &mut F,
) {
    let val = v[last];
    let mut next = last - 1;
    while comp(&val, &v[next]) {
        v[last] = v[next];
        last = next;
        // C++ `--next` past `first` is never dereferenced: the guard element
        // stops the loop first.
        next = next.wrapping_sub(1);
    }
    v[last] = val;
}

/// std::__insertion_sort
fn stdcxx_insertion_sort<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if first == last {
        return;
    }
    for i in (first + 1)..last {
        if comp(&v[i], &v[first]) {
            let val = v[i];
            // _GLIBCXX_MOVE_BACKWARD3(first, i, i + 1)
            v.copy_within(first..i, first + 1);
            v[first] = val;
        } else {
            stdcxx_unguarded_linear_insert(v, i, comp);
        }
    }
}

/// std::__final_insertion_sort
fn stdcxx_final_insertion_sort<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if last - first > STDCXX_S_THRESHOLD {
        stdcxx_insertion_sort(v, first, first + STDCXX_S_THRESHOLD, comp);
        // std::__unguarded_insertion_sort
        for i in (first + STDCXX_S_THRESHOLD)..last {
            stdcxx_unguarded_linear_insert(v, i, comp);
        }
    } else {
        stdcxx_insertion_sort(v, first, last, comp);
    }
}

/// std::__push_heap (offsets relative to `first`)
fn stdcxx_push_heap<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole_index: usize,
    top_index: usize,
    value: T,
    comp: &mut F,
) {
    // C++ `(holeIndex - 1) / 2` on a signed distance: holeIndex == 0 gives 0.
    let mut parent = hole_index.saturating_sub(1) / 2;
    while hole_index > top_index && comp(&v[first + parent], &value) {
        v[first + hole_index] = v[first + parent];
        hole_index = parent;
        parent = hole_index.saturating_sub(1) / 2;
    }
    v[first + hole_index] = value;
}

/// std::__adjust_heap (offsets relative to `first`)
fn stdcxx_adjust_heap<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut hole_index: usize,
    len: usize,
    value: T,
    comp: &mut F,
) {
    let top_index = hole_index;
    let mut second_child = hole_index;
    // len >= 1 here (make_heap / sort_heap guards), so (len - 1) / 2 is safe
    while second_child < (len - 1) / 2 {
        second_child = 2 * (second_child + 1);
        if comp(&v[first + second_child], &v[first + (second_child - 1)]) {
            second_child -= 1;
        }
        v[first + hole_index] = v[first + second_child];
        hole_index = second_child;
    }
    if (len & 1) == 0 && second_child == (len - 2) / 2 {
        second_child = 2 * (second_child + 1);
        v[first + hole_index] = v[first + (second_child - 1)];
        hole_index = second_child - 1;
    }
    stdcxx_push_heap(v, first, hole_index, top_index, value, comp);
}

/// std::__make_heap
fn stdcxx_make_heap<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    comp: &mut F,
) {
    if last - first < 2 {
        return;
    }
    let len = last - first;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent];
        stdcxx_adjust_heap(v, first, parent, len, value, comp);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

/// std::__sort_heap (with std::__pop_heap(first, last, last))
fn stdcxx_sort_heap<T: Copy, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    comp: &mut F,
) {
    while last - first > 1 {
        last -= 1;
        let value = v[last];
        v[last] = v[first];
        stdcxx_adjust_heap(v, first, 0, last - first, value, comp);
    }
}

#[cfg(test)]
mod tests {
    //! Expected values are the C++ oracle's: HiGHS `Analyse(rows, ptr, signs,
    //! nb, nullptr, data, perm).run(S)` from the fixture build
    //! (target/highs-oracle/build/lib/libhighs.a, g++ 13 / libstdc++), dumped
    //! through the `Symbolic` getters in this exact text format (doubles as
    //! IEEE-754 bit patterns).
    use super::super::super::data_collector::DataCollector;
    use super::*;

    fn analyse_dump(
        rows: &[Int],
        ptr: &[Int],
        signs: &[Int],
        nb: Int,
        perm: &[Int],
    ) -> Vec<String> {
        let d = DataCollector::new();
        let mut s = Symbolic::new();
        let mut an = Analyse::new(rows, ptr, signs, nb, None, &d, perm);
        assert_eq!(an.run(&mut s), kRetOk);
        let j = |v: &[Int]| v.iter().map(|x| format!(" {}", x)).collect::<String>();
        let mut out = vec![
            format!(
                "n {} sn {} nz {} bs {} depth {} mss {}",
                s.size(),
                s.sn(),
                s.nz(),
                s.blockSize(),
                s.depth(),
                s.maxStackSize()
            ),
            format!(
                "fillin {:x} flops {:x} spops {:x} critops {:x} storage {:x}",
                s.fillin().to_bits(),
                s.flops().to_bits(),
                s.spops().to_bits(),
                s.critops().to_bits(),
                s.storage().to_bits()
            ),
            format!("iperm{}", j(s.iperm())),
            format!(
                "ptr{}",
                s.ptr_vec()
                    .iter()
                    .map(|x| format!(" {}", x))
                    .collect::<String>()
            ),
        ];
        let last = *s.ptr_vec().last().unwrap();
        out.push(format!(
            "rows{}",
            (0..last)
                .map(|i| format!(" {}", s.rows(i)))
                .collect::<String>()
        ));
        out.push(format!("snparent{}", j(s.snParent_vec())));
        out.push(format!("snstart{}", j(s.snStart_vec())));
        out.push(format!("sign{}", j(s.pivotSign())));
        out.push(format!(
            "relcols{}",
            (0..rows.len() as Int)
                .map(|i| format!(" {}", s.relindCols(i)))
                .collect::<String>()
        ));
        let mut rc = String::from("relclique");
        let mut cs = String::from("consec");
        let mut cb = String::from("cbs");
        for sn in 0..s.sn() {
            let ldf = (s.ptr(sn + 1) - s.ptr(sn)) as Int;
            let ldc = ldf - (s.snStart(sn + 1) - s.snStart(sn));
            if s.snParent(sn) != -1 {
                for k in 0..ldc {
                    rc.push_str(&format!(" {}", s.relindClique(sn, k)));
                    cs.push_str(&format!(" {}", s.consecutiveSums(sn, k)));
                }
                rc.push_str(" |");
                cs.push_str(" |");
            }
            let n_blocks = (ldc - 1) / nb + 1;
            for k in 0..=n_blocks {
                cb.push_str(&format!(" {}", s.cliqueBlockStart(sn, k)));
            }
            cb.push_str(" |");
        }
        out.push(rc);
        out.push(cs);
        out.push(cb);
        out
    }

    #[test]
    fn run_2x2_quasidefinite() {
        // [2 1; 1 -3] lower CSC, signs (+, -), identity ordering, nb = 128.
        let expected = [
            "n 2 sn 1 nz 3 bs 128 depth 1 mss 0",
            "fillin 3ff0000000000000 flops 3ff0000000000000 spops 0 critops 3ff0000000000000 storage 4040000000000000",
            "iperm 0 1",
            "ptr 0 2",
            "rows 0 1",
            "snparent -1",
            "snstart 0 2",
            "sign 1 -1",
            "relcols 0 1 1",
            "relclique",
            "consec",
            "cbs 0 0 |",
        ];
        let got = analyse_dump(&[0, 1, 1], &[0, 2, 3], &[1, -1], 128, &[0, 1]);
        assert_eq!(got, expected);
    }

    #[test]
    fn run_arrow_20_leaves_tied_children() {
        // Arrow: 20 leaves joined to one root; the 20 children tie in reorderChildren,
        // so the order is libstdc++ std::sort's (a stable sort gives iperm 19 0 1 ...).
        let expected = [
            "n 21 sn 20 nz 41 bs 128 depth 2 mss 19",
            "fillin 3ff0000000000000 flops 4034000000000000 spops 4033000000000000 critops 4000000000000000 storage 407e800000000000",
            "iperm 19 10 18 17 16 15 14 13 12 11 0 9 8 7 6 5 4 3 2 1 20",
            "ptr 0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40",
            "rows 0 20 1 20 2 20 3 20 4 20 5 20 6 20 7 20 8 20 9 20 10 20 11 20 12 20 13 20 14 20 15 20 16 20 17 20 18 20 19 20",
            "snparent 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 -1",
            "snstart 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21",
            "sign 1 1 1 -1 1 1 -1 1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 1 -1",
            "relcols 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1",
            "relclique 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |",
            "consec 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |",
            "cbs 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 1 | 0 0 |",
        ];
        let got = analyse_dump(
            &[
                0, 20, 1, 20, 2, 20, 3, 20, 4, 20, 5, 20, 6, 20, 7, 20, 8, 20, 9, 20, 10, 20, 11,
                20, 12, 20, 13, 20, 14, 20, 15, 20, 16, 20, 17, 20, 18, 20, 19, 20, 20,
            ],
            &[
                0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 41,
            ],
            &[
                1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1,
            ],
            128,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ],
        );
        assert_eq!(got, expected);
    }

    #[test]
    fn run_random_nb2() {
        // Random 23x23 pattern, random ordering, nb = 2 (multi-block cliques).
        let expected = [
            "n 23 sn 5 nz 105 bs 2 depth 2 mss 7",
            "fillin 4001df51b3bea367 flops 407ae00000000000 spops 4018000000000000 critops 4075300000000000 storage 408e800000000000",
            "iperm 20 2 15 11 10 8 5 13 6 22 0 19 18 17 1 4 14 7 12 21 16 9 3",
            "ptr 0 1 2 11 19 26",
            "rows 0 1 2 3 4 5 6 7 11 13 15 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22",
            "snparent -1 -1 3 -1 -1",
            "snstart 0 1 2 8 16 23",
            "sign -1 1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 -1 -1 -1 1 -1 1 -1 1 -1",
            "relcols 0 0 0 5 6 8 1 5 7 2 5 3 4 5 6 4 8 5 0 1 2 1 4 2 3 6 3 7 4 5 5 6 7 0 3 6 1 3 6 2 3 5 3 4 4 5 6",
            "relclique 3 5 7 |",
            "consec 1 1 1 |",
            "cbs 0 0 | 0 0 | 0 6 7 | 0 0 | 0 0 |",
        ];
        let got = analyse_dump(
            &[
                0, 11, 1, 2, 3, 17, 2, 3, 8, 3, 4, 6, 4, 5, 16, 5, 21, 6, 8, 17, 7, 18, 22, 8, 9,
                13, 20, 10, 11, 12, 13, 20, 12, 19, 13, 14, 15, 17, 16, 17, 22, 18, 21, 19, 20, 21,
                22,
            ],
            &[
                0, 2, 6, 9, 12, 15, 17, 20, 23, 24, 27, 28, 32, 34, 35, 36, 38, 39, 41, 43, 44, 45,
                46, 47,
            ],
            &[
                -1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1,
            ],
            2,
            &[
                12, 13, 10, 6, 8, 15, 5, 22, 20, 1, 11, 21, 14, 4, 17, 0, 3, 18, 7, 16, 2, 19, 9,
            ],
        );
        assert_eq!(got, expected);
    }

    #[test]
    fn libstdcxx_sort_tie_order() {
        // g++ 13: std::sort of {(i, (i * 7) % 3)} with `a.second > b.second`.
        let expect: [(usize, &[Int]); 3] = [
            (
                17,
                &[14, 2, 11, 5, 8, 16, 13, 10, 7, 4, 1, 9, 6, 12, 3, 15, 0],
            ),
            (
                40,
                &[
                    29, 20, 17, 23, 14, 26, 11, 8, 32, 5, 35, 2, 38, 31, 1, 34, 28, 25, 37, 22, 16,
                    13, 4, 10, 19, 7, 30, 39, 0, 3, 36, 6, 33, 9, 12, 27, 15, 24, 18, 21,
                ],
            ),
            (
                100,
                &[
                    35, 20, 77, 23, 26, 71, 29, 68, 32, 65, 74, 62, 38, 59, 41, 56, 44, 53, 47, 50,
                    89, 98, 2, 95, 5, 92, 8, 80, 11, 86, 14, 17, 83, 79, 70, 52, 82, 97, 55, 76,
                    58, 94, 1, 61, 91, 64, 73, 85, 67, 88, 34, 19, 22, 16, 49, 25, 28, 13, 31, 10,
                    37, 7, 40, 46, 4, 43, 99, 81, 18, 0, 3, 96, 84, 15, 87, 6, 12, 75, 90, 9, 93,
                    36, 51, 48, 54, 45, 57, 42, 60, 39, 63, 21, 66, 33, 69, 30, 72, 27, 24, 78,
                ],
            ),
        ];
        for (n, want) in expect {
            let mut v: Vec<(Int, f64)> = (0..n as Int).map(|i| (i, ((i * 7) % 3) as f64)).collect();
            libstdcxx_sort(&mut v, |a: &(Int, f64), b: &(Int, f64)| a.1 > b.1);
            let got: Vec<Int> = v.iter().map(|p| p.0).collect();
            assert_eq!(got, want, "n = {}", n);
        }
    }
}
