//! HiGHS source: ipm/basiclu/lu_solve_triangular.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_solve_triangular.c::lu_solve_triangular (lines 35-140)
///
/// `end` and `pivot` may be `NULL` in C.
pub fn lu_solve_triangular(
    nz_symb: LuInt,
    pattern_symb: &[LuInt],
    begin: &[LuInt],
    end: Option<&[LuInt]>,
    index: &[LuInt],
    value: &[f64],
    pivot: Option<&[f64]>,
    droptol: f64,
    lhs: &mut [f64],
    pattern: &mut [LuInt],
    flops: &mut LuInt,
) -> LuInt {
    let mut nz: LuInt = 0;
    let mut flop_count: LuInt = 0;

    // rusolve: bounds checks removed in the hot loops. The pattern,
    // column and row indices come from the BASICLU factors, which keep them inside `lhs`,
    // `begin`, `end`, `pivot`, `index` and `value` (checked below by debug_assert in debug
    // builds); the arithmetic and its order are unchanged, so results are bit-identical.
    debug_assert!(pattern_symb.iter().all(|&i| (i as usize) < lhs.len()));
    macro_rules! ix {
        ($a:expr, $i:expr) => {{
            debug_assert!(($i as usize) < $a.len());
            // SAFETY: see the comment above; debug builds assert the bound.
            unsafe { *$a.get_unchecked($i as usize) }
        }};
    }
    macro_rules! lhs_at {
        ($i:expr) => {{
            debug_assert!(($i as usize) < lhs.len());
            // SAFETY: see the comment above; debug builds assert the bound.
            unsafe { lhs.get_unchecked_mut($i as usize) }
        }};
    }

    // The four C branches, in the C order: (pivot && end), (pivot), (end), ().
    for n in 0..nz_symb {
        let ipivot = ix!(pattern_symb, n);
        if *lhs_at!(ipivot) == 0.0 {
            continue;
        }
        let x;
        if let Some(pivot) = pivot {
            *lhs_at!(ipivot) /= ix!(pivot, ipivot);
            x = *lhs_at!(ipivot);
            flop_count += 1;
        } else {
            x = *lhs_at!(ipivot);
        }
        let mut pos = ix!(begin, ipivot);
        match end {
            Some(end) => {
                let e = ix!(end, ipivot);
                while pos < e {
                    let i = ix!(index, pos);
                    *lhs_at!(i) -= x * ix!(value, pos);
                    flop_count += 1;
                    pos += 1;
                }
            }
            None => loop {
                let i = ix!(index, pos);
                if i < 0 {
                    break;
                }
                *lhs_at!(i) -= x * ix!(value, pos);
                flop_count += 1;
                pos += 1;
            },
        }
        if x.abs() > droptol {
            pattern[nz as usize] = ipivot;
            nz += 1;
        } else {
            *lhs_at!(ipivot) = 0.0;
        }
    }

    *flops += flop_count;
    nz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_variant_without_pivot_keeps_all_entries() {
        // column 0 holds rows 1 and 2 with values 2 and 3; lhs[0] = 1
        let begin: [LuInt; 3] = [0, 2, 2];
        let end: [LuInt; 3] = [2, 2, 2];
        let index: [LuInt; 2] = [1, 2];
        let value = [2.0, 3.0];
        let pattern_symb = [0];
        let mut lhs = [1.0, 0.0, 0.0];
        let mut pattern = [0; 3];
        let mut flops = 0;
        let nz = lu_solve_triangular(
            1,
            &pattern_symb,
            &begin,
            Some(&end),
            &index,
            &value,
            None,
            0.0,
            &mut lhs,
            &mut pattern,
            &mut flops,
        );
        assert_eq!(nz, 1);
        assert_eq!(pattern[0], 0);
        assert_eq!(lhs, [1.0, -2.0, -3.0]);
        assert_eq!(flops, 2);
    }

    #[test]
    fn pivot_division_then_droptol_drops_small_multiplier() {
        // pivot 2 gives x = 0.5, which is <= droptol 0.6: lhs[0] is zeroed,
        // the pattern is empty, and the update of rows 1 and 2 still happens.
        let begin: [LuInt; 3] = [0, 2, 2];
        let end: [LuInt; 3] = [2, 2, 2];
        let index: [LuInt; 2] = [1, 2];
        let value = [2.0, 3.0];
        let pivot = [2.0, 1.0, 1.0];
        let pattern_symb = [0];
        let mut lhs = [1.0, 0.0, 0.0];
        let mut pattern = [0; 3];
        let mut flops = 0;
        let nz = lu_solve_triangular(
            1,
            &pattern_symb,
            &begin,
            Some(&end),
            &index,
            &value,
            Some(&pivot),
            0.6,
            &mut lhs,
            &mut pattern,
            &mut flops,
        );
        assert_eq!(nz, 0);
        assert_eq!(lhs, [0.0, -1.0, -1.5]);
        assert_eq!(flops, 3);
    }
}
