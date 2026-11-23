//! HiGHS source: ipm/basiclu/lu_list.h (static inline doubly linked lists)
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Shared by both BASICLU factorization paths.
//! Elements `0..nelem` are list members; `nelem + list` is the head of
//! `list`. C `assert`s become `debug_assert!`.

use super::*;

/// HiGHS: basiclu/lu_list.h::lu_list_init (line 56)
pub fn lu_list_init(
    flink: &mut [LuInt],
    blink: &mut [LuInt],
    nelem: LuInt,
    nlist: LuInt,
    min_list: Option<&mut LuInt>,
) {
    for i in 0..nelem + nlist {
        flink[i as usize] = i;
        blink[i as usize] = i;
    }
    if let Some(min_list) = min_list {
        *min_list = 1.max(nlist);
    }
}

/// HiGHS: basiclu/lu_list.h::lu_list_add (line 73)
pub fn lu_list_add(
    elem: LuInt,
    list: LuInt,
    flink: &mut [LuInt],
    blink: &mut [LuInt],
    nelem: LuInt,
    min_list: Option<&mut LuInt>,
) {
    let e = elem as usize;
    let h = (nelem + list) as usize;
    debug_assert!(flink[e] == elem);
    debug_assert!(blink[e] == elem);
    // append elem to the end of list
    let temp = blink[h];
    blink[h] = elem;
    blink[e] = temp;
    flink[temp as usize] = elem;
    flink[e] = nelem + list;
    if let Some(min_list) = min_list {
        if list > 0 && list < *min_list {
            *min_list = list;
        }
    }
}

/// HiGHS: basiclu/lu_list.h::lu_list_remove (line 98)
pub fn lu_list_remove(flink: &mut [LuInt], blink: &mut [LuInt], elem: LuInt) {
    let e = elem as usize;
    flink[blink[e] as usize] = flink[e];
    blink[flink[e] as usize] = blink[e];
    flink[e] = elem;
    blink[e] = elem;
}

/// HiGHS: basiclu/lu_list.h::lu_list_move (line 114)
pub fn lu_list_move(
    elem: LuInt,
    list: LuInt,
    flink: &mut [LuInt],
    blink: &mut [LuInt],
    nelem: LuInt,
    min_list: Option<&mut LuInt>,
) {
    lu_list_remove(flink, blink, elem);
    lu_list_add(elem, list, flink, blink, nelem, min_list);
}

/// HiGHS: basiclu/lu_list.h::lu_list_swap (line 131)
pub fn lu_list_swap(flink: &mut [LuInt], blink: &mut [LuInt], e1: LuInt, e2: LuInt) {
    let (u1, u2) = (e1 as usize, e2 as usize);
    let e1next = flink[u1];
    let e2next = flink[u2];
    let e1prev = blink[u1];
    let e2prev = blink[u2];

    debug_assert!(e1next != e1); // must be in a list
    debug_assert!(e2next != e2);

    if e1next == e2 {
        flink[u2] = e1;
        blink[u1] = e2;
        flink[e1prev as usize] = e2;
        blink[u2] = e1prev;
        flink[u1] = e2next;
        blink[e2next as usize] = e1;
    } else if e2next == e1 {
        flink[u1] = e2;
        blink[u2] = e1;
        flink[u2] = e1next;
        blink[e1next as usize] = e2;
        flink[e2prev as usize] = e1;
        blink[u1] = e2prev;
    } else {
        flink[u2] = e1next;
        blink[e1next as usize] = e2;
        flink[e2prev as usize] = e1;
        blink[u1] = e2prev;
        flink[e1prev as usize] = e2;
        blink[u2] = e1prev;
        flink[u1] = e2next;
        blink[e2next as usize] = e1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn list_add_remove_swap_roundtrip() {
        // 4 elements, 2 lists; heads at 4 and 5.
        let mut f = [0; 6];
        let mut b = [0; 6];
        let mut min_list = 0;
        lu_list_init(&mut f, &mut b, 4, 2, Some(&mut min_list));
        assert_eq!(min_list, 2);
        lu_list_add(0, 1, &mut f, &mut b, 4, Some(&mut min_list));
        lu_list_add(2, 1, &mut f, &mut b, 4, Some(&mut min_list));
        lu_list_add(3, 1, &mut f, &mut b, 4, Some(&mut min_list));
        assert_eq!(min_list, 1);
        // list 1: 5 -> 0 -> 2 -> 3 -> 5
        assert_eq!((f[5], f[0], f[2], f[3]), (0, 2, 3, 5));
        lu_list_swap(&mut f, &mut b, 0, 3);
        assert_eq!((f[5], f[3], f[2], f[0]), (3, 2, 0, 5));
        lu_list_remove(&mut f, &mut b, 2);
        assert_eq!((f[5], f[3], f[0]), (3, 0, 5));
        assert_eq!((b[5], b[0], b[3]), (0, 3, 5));
        assert_eq!((f[2], b[2]), (2, 2));
    }
}
