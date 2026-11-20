//! Index-driven work vector for hyper-sparse simplex passes.
//!
//! A full-length value array plus the list of positions that may hold a
//! nonzero — the shape of HiGHS's `HVector` (`util/HVectorBase.h`: `count`,
//! `index`, `array`; MIT, read as a reference text and reimplemented, not
//! translated) and of the "index-driven" vectors of Hall & McKinnon,
//! *Hyper-sparsity in the revised simplex method and how to exploit it*,
//! Comput. Optim. Appl. 32:259–283 (2005). Every pass that consumes one
//! (solve entry, result scatter, x_B / DSE / δ updates, the clear itself)
//! costs O(pattern) instead of O(n) while the pattern is valid.
//!
//! Differences from `HVector`, both for bit-exactness against the dense
//! loops this type replaces:
//!
//! * No tiny-value dropping (`kHighsTiny`) and no `kHighsZero` placeholder:
//!   values are exactly what the dense code computes. Duplicate-free
//!   patterns come from a `mark` array instead (an entry that cancels to
//!   0.0 stays listed — the pattern is a superset of the nonzeros).
//! * Explicit **dense mode**: once the pattern would exceed
//!   [`HVec::cap`] (or a producer cannot know its pattern) the vector
//!   forgets it and every consumer falls back to its O(n) loop. HiGHS makes
//!   the same switch with its density thresholds (`kHyperCancel` /
//!   `kHyperResult`, `util/HFactorConst.h`); ours is `n / 10`.
//!
//! Invariant (sparse mode): `val[i] != 0.0` ⇒ `mark[i]` ⇒ `i ∈ idx`, and
//! `idx` has no duplicates. In dense mode `idx` is empty and `mark` is all
//! false.

#[derive(Debug, Clone, Default)]
pub(crate) struct HVec {
    val: Vec<f64>,
    idx: Vec<u32>,
    mark: Vec<bool>,
    dense: bool,
    /// `idx` is known to be ascending (set by [`HVec::sort`], cleared by any
    /// append that breaks it).
    sorted: bool,
    cap: usize,
}

impl HVec {
    pub(crate) fn new(n: usize) -> Self {
        HVec {
            val: vec![0.0; n],
            idx: Vec::new(),
            mark: vec![false; n],
            dense: false,
            sorted: true,
            cap: (n / 10).max(1),
        }
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.val.len()
    }

    /// Resize to `n` (all zero, sparse, empty pattern) when the length
    /// differs; a no-op otherwise.
    pub(crate) fn ensure_len(&mut self, n: usize) {
        if self.val.len() != n {
            *self = HVec::new(n);
        }
    }

    #[inline]
    pub(crate) fn is_dense(&self) -> bool {
        self.dense
    }

    /// The value array (always the true vector, in either mode).
    #[inline]
    pub(crate) fn vals(&self) -> &[f64] {
        &self.val
    }

    /// Pattern (superset of the nonzero positions), `None` in dense mode.
    #[inline]
    pub(crate) fn pattern(&self) -> Option<&[u32]> {
        if self.dense { None } else { Some(&self.idx) }
    }

    /// Number of listed positions (sparse mode) or `n` (dense).
    #[inline]
    pub(crate) fn count(&self) -> usize {
        if self.dense {
            self.val.len()
        } else {
            self.idx.len()
        }
    }

    /// Zero the vector: O(pattern) in sparse mode, O(n) in dense mode.
    /// Returns it to sparse mode with an empty pattern.
    pub(crate) fn clear(&mut self) {
        if self.dense {
            self.val.fill(0.0);
            self.dense = false;
        } else {
            for &i in &self.idx {
                self.val[i as usize] = 0.0;
                self.mark[i as usize] = false;
            }
        }
        self.idx.clear();
        self.sorted = true;
    }

    /// Forget the pattern (values untouched): every consumer runs dense.
    pub(crate) fn make_dense(&mut self) {
        if !self.dense {
            for &i in &self.idx {
                self.mark[i as usize] = false;
            }
            self.idx.clear();
            self.dense = true;
        }
    }

    /// Dense write access: drops the pattern first (the writer may touch
    /// any entry).
    #[inline]
    pub(crate) fn dense_mut(&mut self) -> &mut [f64] {
        self.make_dense();
        &mut self.val
    }

    /// Record that position `i` may become nonzero. Switches to dense mode
    /// when the pattern outgrows `cap`.
    #[inline]
    pub(crate) fn touch(&mut self, i: usize) {
        if !self.dense && !self.mark[i] {
            if self.idx.len() >= self.cap {
                self.make_dense();
                return;
            }
            self.mark[i] = true;
            if self.sorted && self.idx.last().is_some_and(|&l| l as usize > i) {
                self.sorted = false;
            }
            self.idx.push(i as u32);
        }
    }

    /// `val[i] = v` with pattern bookkeeping.
    #[inline]
    pub(crate) fn set(&mut self, i: usize, v: f64) {
        self.touch(i);
        self.val[i] = v;
    }

    /// Raw value write for a position the caller has already listed (or for
    /// any position in dense mode). Debug-checked.
    #[inline]
    pub(crate) fn put_listed(&mut self, i: usize, v: f64) {
        debug_assert!(self.dense || self.mark[i] || v == 0.0);
        self.val[i] = v;
    }

    /// `val[i] = v`, listing `i` only when `v != 0.0`. An unlisted position
    /// already holds +0.0, so skipping a zero write there changes at most
    /// the sign of a zero — this is what keeps the "unlisted ⇒ +0.0"
    /// invariant that [`HVec::clear`] relies on.
    #[inline]
    pub(crate) fn assign(&mut self, i: usize, v: f64) {
        if v != 0.0 {
            self.touch(i);
            self.val[i] = v;
        } else if self.dense || self.mark[i] {
            self.val[i] = v;
        }
    }

    /// [`HVec::touch`] keeping an ascending pattern ascending (sorted
    /// insertion, O(count)): for the BTRAN eta walk, whose dot products
    /// binary-search the pattern against the eta's ascending index list.
    #[inline]
    pub(crate) fn touch_sorted(&mut self, i: usize) {
        if !self.dense && !self.mark[i] {
            if self.idx.len() >= self.cap {
                self.make_dense();
                return;
            }
            self.mark[i] = true;
            if self.sorted {
                let k = self.idx.partition_point(|&j| (j as usize) < i);
                self.idx.insert(k, i as u32);
            } else {
                self.idx.push(i as u32);
            }
        }
    }

    /// Ascending `(i, v)` list of the nonzeros (sorts the pattern; O(n) scan
    /// in dense mode) — the sparse column shape `basis_nz` stores.
    pub(crate) fn nonzeros_sorted(&mut self) -> Vec<(usize, f64)> {
        self.sort();
        match self.pattern() {
            Some(p) => p
                .iter()
                .filter_map(|&i| {
                    let v = self.val[i as usize];
                    (v != 0.0).then_some((i as usize, v))
                })
                .collect(),
            None => self
                .val
                .iter()
                .enumerate()
                .filter_map(|(i, &v)| (v != 0.0).then_some((i, v)))
                .collect(),
        }
    }

    /// Mutable value of position `i`, listing it first (value semantics as
    /// for a plain `&mut f64`: the caller may write anything).
    #[inline]
    pub(crate) fn entry(&mut self, i: usize) -> &mut f64 {
        self.touch(i);
        &mut self.val[i]
    }

    /// Sort the pattern ascending (no-op when already sorted / dense):
    /// consumers whose fp accumulation order is "ascending index" need it to
    /// stay bit-identical to their dense loops.
    pub(crate) fn sort(&mut self) {
        if !self.dense && !self.sorted {
            self.idx.sort_unstable();
            self.sorted = true;
        }
    }

    /// Copy `src` into `self` (same length), O(pattern) when both are sparse.
    pub(crate) fn copy_from(&mut self, src: &HVec) {
        debug_assert_eq!(self.len(), src.len());
        self.clear();
        match src.pattern() {
            Some(p) if p.len() <= self.cap => {
                for &i in p {
                    let i = i as usize;
                    self.mark[i] = true;
                    self.val[i] = src.val[i];
                }
                self.idx.extend_from_slice(p);
                self.sorted = src.sorted;
            }
            _ => {
                self.val.copy_from_slice(&src.val);
                self.dense = true;
            }
        }
    }

    /// Load a dense slice (dense mode).
    pub(crate) fn load_dense(&mut self, src: &[f64]) {
        self.make_dense();
        self.val.copy_from_slice(src);
    }

    /// Recover a pattern for a dense-mode vector by scanning (O(n)); stays
    /// dense when the nonzero count exceeds `cap`. Test-only today.
    #[cfg(test)]
    pub(crate) fn reindex(&mut self) {
        if !self.dense {
            return;
        }
        let mut cnt = 0usize;
        for &v in &self.val {
            cnt += (v != 0.0) as usize;
        }
        if cnt > self.cap {
            return;
        }
        self.dense = false;
        self.idx.clear();
        for (i, &v) in self.val.iter().enumerate() {
            if v != 0.0 {
                self.mark[i] = true;
                self.idx.push(i as u32);
            }
        }
        self.sorted = true;
    }

    /// Max-abs over the vector, O(pattern) when sparse.
    pub(crate) fn max_abs(&self) -> f64 {
        match self.pattern() {
            Some(p) => p
                .iter()
                .fold(0.0f64, |a, &i| a.max(self.val[i as usize].abs())),
            None => self.val.iter().fold(0.0f64, |a, &v| a.max(v.abs())),
        }
    }

    /// Number of nonzero values, O(pattern) when sparse.
    pub(crate) fn nnz(&self) -> usize {
        match self.pattern() {
            Some(p) => p.iter().filter(|&&i| self.val[i as usize] != 0.0).count(),
            None => self.val.iter().filter(|v| **v != 0.0).count(),
        }
    }

    /// Pattern capacity before the vector turns dense.
    #[inline]
    pub(crate) fn cap(&self) -> usize {
        self.cap
    }
}

impl std::ops::Index<usize> for HVec {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        &self.val[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pattern_bookkeeping_and_dense_switch() {
        let mut h = HVec::new(40); // cap 4
        h.set(7, 1.0);
        h.set(3, 2.0);
        h.set(7, 3.0);
        assert_eq!(h.pattern().unwrap(), &[7, 3]);
        h.sort();
        assert_eq!(h.pattern().unwrap(), &[3, 7]);
        assert_eq!(h.nnz(), 2);
        assert_eq!(h.max_abs(), 3.0);
        h.set(1, 1.0);
        h.set(2, 1.0);
        assert!(!h.is_dense());
        h.set(30, 5.0); // 5th entry > cap ⇒ dense, value still written
        assert!(h.is_dense());
        assert_eq!(h[30], 5.0);
        assert_eq!(h.nnz(), 5);
        h.clear();
        assert!(!h.is_dense());
        assert!(h.vals().iter().all(|&v| v == 0.0));
        h.set(5, 1.0);
        assert_eq!(h.pattern().unwrap(), &[5]);
        let mut g = HVec::new(40);
        g.copy_from(&h);
        assert_eq!(g.pattern().unwrap(), &[5]);
        g.dense_mut()[9] = 2.0;
        g.reindex();
        assert_eq!(g.pattern().unwrap(), &[5, 9]);
    }
}
