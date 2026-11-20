//! HiGHS source: ipm/basiclu/lu_dfs.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_dfs.c::lu_dfs (lines 42-52)
///
/// `end` is `NULL` for the `dfs` variant. `pstack` is the `lu_int` stack
/// punned over `work1` ([`pun_set_int`]/[`pun_get_int`]).
pub fn lu_dfs(
    i: LuInt,
    begin: &[LuInt],
    end: Option<&[LuInt]>,
    index: &[LuInt],
    top: LuInt,
    xi: &mut [LuInt],
    pstack: &mut [f64],
    marked: &mut [LuInt],
    M: LuInt,
) -> LuInt {
    if marked[i as usize] == M {
        return top;
    }

    match end {
        Some(end) => dfs_end(i, begin, end, index, top, xi, pstack, marked, M),
        None => dfs(i, begin, index, top, xi, pstack, marked, M),
    }
}

// rusolve: bounds checks removed in the two DFS loops. Node numbers
// and positions come from the BASICLU factors and lie inside `begin`, `end`, `index`, `marked`;
// the stack depth `head` and `top` stay inside `xi` and `pstack` (xi and pstack have length m
// and every node is pushed at most once, guarded by `marked`). Debug builds assert each bound.
// The visiting order and the result are unchanged.
macro_rules! ix {
    ($a:expr, $i:expr) => {{
        debug_assert!(($i as usize) < $a.len());
        // SAFETY: see the comment above; debug builds assert the bound.
        unsafe { *$a.get_unchecked($i as usize) }
    }};
}
macro_rules! put {
    ($a:expr, $i:expr, $v:expr) => {{
        debug_assert!(($i as usize) < $a.len());
        // SAFETY: see the comment above; debug builds assert the bound.
        unsafe { *$a.get_unchecked_mut($i as usize) = $v }
    }};
}

/// HiGHS: basiclu/lu_dfs.c::dfs_end (lines 61-97)
fn dfs_end(
    i: LuInt,
    begin: &[LuInt],
    end: &[LuInt],
    index: &[LuInt],
    top: LuInt,
    xi: &mut [LuInt],
    pstack: &mut [f64],
    marked: &mut [LuInt],
    M: LuInt,
) -> LuInt {
    let mut i = i;
    let mut top = top;
    let mut head: LuInt = 0;
    debug_assert!(marked[i as usize] != M);

    xi[0] = i;
    while head >= 0 {
        i = ix!(xi, head);
        if ix!(marked, i) != M {
            // node i has not been visited
            put!(marked, i, M);
            let b = ix!(begin, i);
            put!(pstack, head, f64::from_bits(b as i64 as u64));
        }
        let mut done: LuInt = 1;
        // continue dfs at node i
        let mut p = ix!(pstack, head).to_bits() as i64 as LuInt;
        let pend = ix!(end, i);
        while p < pend {
            let inext = ix!(index, p);
            if ix!(marked, inext) == M {
                // skip visited node
                p += 1;
                continue;
            }
            put!(pstack, head, f64::from_bits((p + 1) as i64 as u64));
            head += 1;
            put!(xi, head, inext); // start dfs at node inext
            done = 0;
            break;
        }
        if done != 0 {
            // node i has no unvisited neighbours
            head -= 1;
            top -= 1;
            put!(xi, top, i);
        }
    }

    top
}

/// HiGHS: basiclu/lu_dfs.c::dfs (lines 106-141)
fn dfs(
    i: LuInt,
    begin: &[LuInt],
    index: &[LuInt],
    top: LuInt,
    xi: &mut [LuInt],
    pstack: &mut [f64],
    marked: &mut [LuInt],
    M: LuInt,
) -> LuInt {
    let mut i = i;
    let mut top = top;
    let mut head: LuInt = 0;
    debug_assert!(marked[i as usize] != M);

    xi[0] = i;
    while head >= 0 {
        i = ix!(xi, head);
        if ix!(marked, i) != M {
            // node i has not been visited
            put!(marked, i, M);
            let b = ix!(begin, i);
            put!(pstack, head, f64::from_bits(b as i64 as u64));
        }
        let mut done: LuInt = 1;
        // continue dfs at node i; the list is terminated by a negative index
        let mut p = ix!(pstack, head).to_bits() as i64 as LuInt;
        loop {
            let inext = ix!(index, p);
            if inext < 0 {
                break;
            }
            if ix!(marked, inext) == M {
                // skip visited node
                p += 1;
                continue;
            }
            put!(pstack, head, f64::from_bits((p + 1) as i64 as u64));
            head += 1;
            put!(xi, head, inext); // start dfs at node inext
            done = 0;
            break;
        }
        if done != 0 {
            // node i has no unvisited neighbours
            head -= 1;
            top -= 1;
            put!(xi, top, i);
        }
    }

    top
}

#[cfg(test)]
mod tests {
    use super::*;

    // Graph on 4 nodes: 0 -> [1, 2], 1 -> [3], 2 -> [], 3 -> [].
    // Hand-traced DFS from node 0 gives postorder 3, 1, 2, 0 written from
    // the top of xi downwards, so xi[0..4] = [0, 2, 1, 3] and top = 0.

    #[test]
    fn dfs_end_variant_reach_order() {
        let begin: [LuInt; 4] = [0, 2, 3, 3];
        let end: [LuInt; 4] = [2, 3, 3, 3];
        let index: [LuInt; 3] = [1, 2, 3];
        let mut xi = [0; 4];
        let mut pstack = [0.0f64; 4];
        let mut marked = [0; 4];
        let top = lu_dfs(
            0,
            &begin,
            Some(&end),
            &index,
            4,
            &mut xi,
            &mut pstack,
            &mut marked,
            1,
        );
        assert_eq!(top, 0);
        assert_eq!(xi, [0, 2, 1, 3]);
        assert_eq!(marked, [1, 1, 1, 1]);
    }

    #[test]
    fn dfs_sentinel_variant_reach_order() {
        // same graph, neighbour lists terminated by -1
        let begin: [LuInt; 4] = [0, 3, 5, 6];
        let index: [LuInt; 7] = [1, 2, -1, 3, -1, -1, -1];
        let mut xi = [0; 4];
        let mut pstack = [0.0f64; 4];
        let mut marked = [0; 4];
        let top = lu_dfs(
            0,
            &begin,
            None,
            &index,
            4,
            &mut xi,
            &mut pstack,
            &mut marked,
            1,
        );
        assert_eq!(top, 0);
        assert_eq!(xi, [0, 2, 1, 3]);
    }

    #[test]
    fn dfs_marked_start_is_noop() {
        let begin: [LuInt; 4] = [0, 2, 3, 3];
        let end: [LuInt; 4] = [2, 3, 3, 3];
        let index: [LuInt; 3] = [1, 2, 3];
        let mut xi = [0; 4];
        let mut pstack = [0.0f64; 4];
        let mut marked = [5, 0, 0, 0];
        let top = lu_dfs(
            0,
            &begin,
            Some(&end),
            &index,
            4,
            &mut xi,
            &mut pstack,
            &mut marked,
            5,
        );
        assert_eq!(top, 4);
    }
}
