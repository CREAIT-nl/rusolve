//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::*;

#[derive(Clone, Default, Debug)]
pub(crate) struct S6SolPool {
    pub sols: Vec<S6Sol>,
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsLT (lines 6866-6883), via def.h EPSLT (line 183)
#[inline]
fn s6_set_is_lt(a: f64, b: f64) -> bool {
    (a - b) < -S6_EPSILON
}

/// SCIP 11.0.0 (32dd425): set.c::SCIPsetIsGT (lines 6906-6923), via def.h EPSGT (line 185)
#[inline]
fn s6_set_is_gt(a: f64, b: f64) -> bool {
    (a - b) > S6_EPSILON
}

impl S6SolPool {
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.sols.len()
    }

    #[inline]
    pub(crate) fn best(&self) -> Option<&S6Sol> {
        self.sols.first()
    }

    /// SCIP 11.0.0 (32dd425): primal.c::primalSearchSolPos (lines 1205-1250)
    ///
    /// The transformed-space preference loop (primal.c:1243-1248) is omitted
    /// because every pool sol is transformed (no original-space sols exist).
    pub(crate) fn s6_search_sol_pos(&self, obj: f64) -> usize {
        let sols = &self.sols;
        let nsols = sols.len() as i64;

        let mut left: i64 = -1;
        let mut right: i64 = nsols;
        while left < right - 1 {
            // C truncating division: (left+right)/2
            let middle = (left + right) / 2;
            debug_assert!(left < middle && middle < right);
            debug_assert!(0 <= middle && middle < nsols);

            let middleobj = sols[middle as usize].obj;

            if obj < middleobj {
                right = middle;
            } else {
                left = middle;
            }
        }
        debug_assert!(left == right - 1);

        // preference loop omitted (see doc comment)

        right as usize
    }

    /// SCIP 11.0.0 (32dd425): primal.c::primalExistsSol (lines 1290-1368)
    ///
    /// The replace branch (origin original) is omitted because every pool sol is
    /// transformed (no original-space sols exist).
    pub(crate) fn s6_exists_sol(&self, obj: f64, x: &[f64], insertpos: usize) -> bool {
        debug_assert!(insertpos <= self.sols.len());

        // search in the better solutions
        let mut i: i64 = insertpos as i64 - 1;
        while i >= 0 {
            let solobj = self.sols[i as usize].obj;

            if s6_set_is_lt(solobj, obj) {
                break;
            }

            if s6_sols_are_equal(obj, x, solobj, &self.sols[i as usize].x) {
                return true;
            }

            i -= 1;
        }

        // search in the worse solutions
        let mut i: usize = insertpos;
        while i < self.sols.len() {
            let solobj = self.sols[i].obj;

            if s6_set_is_gt(solobj, obj) {
                break;
            }

            if s6_sols_are_equal(obj, x, solobj, &self.sols[i].x) {
                return true;
            }

            i += 1;
        }

        false
    }

    /// SCIP 11.0.0 (32dd425): primal.c::solOfInterest (lines 1470-1497) +
    /// primal.c::primalAddSol (lines 859-1100, non-replace path)
    ///
    /// misc_improvingsols is FALSE (set.c:323), so the cutoff test is not applied.
    /// The stats, events, upper bound and visualisation updates in primalAddSol are
    /// not ported here (caller side, D4).
    pub(crate) fn insert(
        &mut self,
        obj: f64,
        x: Vec<f64>,
        origin: S6Origin,
        nodenum: i64,
        next_index: &mut i64,
    ) -> Option<(usize, i64)> {
        // solOfInterest :1473-1479: infinite objective is rejected
        if s6_is_inf(obj) {
            return None;
        }

        // solOfInterest :1486
        let pos = self.s6_search_sol_pos(obj);

        // solOfInterest :1492-1493: short-circuit as in C
        if pos < S6_MAXSOL && !self.s6_exists_sol(obj, &x, pos) {
            // primalAddSol :1021-1024: if pool is full, drop the last sol
            if self.sols.len() == S6_MAXSOL {
                self.sols.pop();
            }

            // primalAddSol :1032-1037: shift and insert (Vec::insert does the shift)
            let index = *next_index;
            *next_index += 1;

            let sol = S6Sol {
                index,
                obj,
                x,
                origin,
                nodenum,
            };
            self.sols.insert(pos, sol);

            Some((pos, index))
        } else {
            None
        }
    }
}

/// SCIP 11.0.0 (32dd425): sol.c::SCIPsolsAreEqual (lines 3385-3452)
///
/// Transformed-space only. The exact-mode branch (set->exact_enable) is not ported.
/// NOTE: the C does NOT return false for infinite objectives of the same sign
/// (infobjs skips the objective test, and the value loop runs); ported exactly.
pub(crate) fn s6_sols_are_equal(obj1: f64, x1: &[f64], obj2: f64, x2: &[f64]) -> bool {
    // sol.c:3427-3428: infobjs
    let infobjs = (s6_is_inf(obj1) && s6_is_inf(obj2)) || (s6_is_inf(-obj1) && s6_is_inf(-obj2));

    // sol.c:3429-3430
    if !infobjs && !s6_is_eq(obj1, obj2) {
        return false;
    }

    // sol.c:3441-3450: compare each variable value (prob->nvars, same for both sols)
    debug_assert_eq!(x1.len(), x2.len());
    for v in 0..x1.len() {
        if !s6_is_eq(x1[v], x2[v]) {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn s6_pool_equal_obj_insert_after_existing() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        let r1 = pool.insert(5.0, vec![1.0, 0.0], S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(r1, Some((0, 0)));

        let r2 = pool.insert(5.0, vec![0.0, 1.0], S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(r2, Some((1, 1)));

        assert_eq!(pool.s6_search_sol_pos(5.0), 2);
        assert_eq!(pool.sols[0].index, 0);
        assert_eq!(pool.sols[1].index, 1);
    }

    #[test]
    fn s6_pool_duplicate_rejected_index_not_consumed() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        pool.insert(5.0, vec![1.0, 0.0], S6Origin::Highs(0), 0, &mut next_index);
        pool.insert(5.0, vec![0.0, 1.0], S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(next_index, 2);

        let r = pool.insert(5.0, vec![1.0, 0.0], S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(r, None);
        assert_eq!(next_index, 2);
    }

    #[test]
    fn s6_pool_full_drops_worst() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        for i in 0..100 {
            let obj = i as f64;
            let r = pool.insert(obj, vec![obj], S6Origin::Highs(0), 0, &mut next_index);
            assert!(r.is_some(), "insert obj {} failed", obj);
        }
        assert_eq!(pool.len(), 100);

        // obj 1000.0: pos = 100, not < S6_MAXSOL
        let r = pool.insert(1000.0, vec![1000.0], S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(r, None);

        // obj -1.0: pos = 0, pool is full so last sol (obj 99.0) is dropped
        let r = pool.insert(-1.0, vec![-1.0], S6Origin::Highs(0), 0, &mut next_index);
        assert!(r.is_some());
        assert_eq!(r.unwrap().0, 0);
        assert_eq!(pool.len(), 100);
        assert_eq!(pool.best().unwrap().obj, -1.0);

        // obj 99.0 was dropped
        let has_99 = pool.sols.iter().any(|s| s.obj == 99.0);
        assert!(!has_99);
    }

    #[test]
    fn s6_pool_eps_equal_objectives() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        let obj1 = 1.0;
        let obj2 = 1.0 + 1e-10;
        let x = vec![1.0, 0.0];

        assert!(s6_sols_are_equal(obj1, &x, obj2, &x));

        pool.insert(obj1, x.clone(), S6Origin::Highs(0), 0, &mut next_index);
        let r = pool.insert(obj2, x.clone(), S6Origin::Highs(0), 0, &mut next_index);
        assert_eq!(r, None);
    }

    #[test]
    fn s6_pool_infinite_objective_rejected() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        let r = pool.insert(
            S6_INFINITY,
            vec![1.0],
            S6Origin::Highs(0),
            0,
            &mut next_index,
        );
        assert_eq!(r, None);
        assert_eq!(next_index, 0);

        // C behaviour: infobjs skips the objective test, value loop runs
        assert!(s6_sols_are_equal(S6_INFINITY, &[1.0], S6_INFINITY, &[1.0]));
        assert!(!s6_sols_are_equal(
            S6_INFINITY,
            &[1.0],
            -S6_INFINITY,
            &[1.0]
        ));
    }

    #[test]
    fn s6_pool_indices_monotone() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        let objs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mut indices = Vec::new();
        for &obj in &objs {
            let r = pool.insert(obj, vec![obj], S6Origin::Highs(0), 0, &mut next_index);
            assert!(r.is_some());
            indices.push(r.unwrap().1);
        }
        assert_eq!(indices, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn s6_pool_search_pos_binary() {
        let mut pool = S6SolPool::default();
        let mut next_index: i64 = 0;

        // insert objs [1.0, 2.0, 2.0, 3.0] in that order
        pool.insert(1.0, vec![1.0], S6Origin::Highs(0), 0, &mut next_index);
        pool.insert(2.0, vec![2.0], S6Origin::Highs(0), 0, &mut next_index);
        pool.insert(2.0, vec![2.5], S6Origin::Highs(0), 0, &mut next_index);
        pool.insert(3.0, vec![3.0], S6Origin::Highs(0), 0, &mut next_index);

        assert_eq!(pool.s6_search_sol_pos(2.0), 3);
        assert_eq!(pool.s6_search_sol_pos(0.5), 0);
        assert_eq!(pool.s6_search_sol_pos(4.0), 4);
    }
}
