//! HiGHS source: ipm/basiclu/lu_solve_symbolic.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU, from the C source header)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_dfs::lu_dfs;
use super::lu_internal::Lu;
use super::*;

/// HiGHS: basiclu/lu_solve_symbolic.c::lu_solve_symbolic (lines 27-48)
pub fn lu_solve_symbolic(
    m: LuInt,
    begin: &[LuInt],
    end: Option<&[LuInt]>,
    index: &[LuInt],
    nrhs: LuInt,
    irhs: &[LuInt],
    ilhs: &mut [LuInt],
    pstack: &mut [f64],
    marked: &mut [LuInt],
    M: LuInt,
) -> LuInt {
    let mut top = m;

    for n in 0..nrhs {
        let i = irhs[n as usize];
        if marked[i as usize] != M {
            top = lu_dfs(i, begin, end, index, top, ilhs, pstack, marked, M);
        }
    }

    top
}
