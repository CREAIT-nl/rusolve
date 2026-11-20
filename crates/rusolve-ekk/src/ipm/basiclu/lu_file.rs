//! HiGHS source: ipm/basiclu/lu_file.c
//!
//! Copyright (C) 2016-2018 ERGO-Code (BASICLU, lu_file.c)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_list::lu_list_move;
use super::*;

/// HiGHS: basiclu/lu_file.c::lu_file_empty (lines 49-66)
pub fn lu_file_empty(
    nlines: LuInt,
    begin: &mut [LuInt],
    end: &mut [LuInt],
    next: &mut [LuInt],
    prev: &mut [LuInt],
    fmem: LuInt,
) {
    let n = nlines as usize;
    begin[n] = 0;
    end[n] = fmem;
    for i in 0..n {
        begin[i] = 0;
        end[i] = 0;
    }
    for i in 0..n {
        next[i] = i as LuInt + 1;
        prev[i + 1] = i as LuInt;
    }
    next[n] = 0;
    prev[0] = nlines;
}

/// HiGHS: basiclu/lu_file.c::lu_file_reappend (lines 76-100)
#[allow(clippy::too_many_arguments)] // C signature (BASICLU file API)
pub fn lu_file_reappend(
    line: LuInt,
    nlines: LuInt,
    begin: &mut [LuInt],
    end: &mut [LuInt],
    next: &mut [LuInt],
    prev: &mut [LuInt],
    index: &mut [LuInt],
    value: &mut [f64],
    extra_space: LuInt,
) {
    let l = line as usize;
    let n = nlines as usize;
    let fmem = end[n];
    let mut used = begin[n];
    let ibeg = begin[l]; // old beginning of line
    let iend = end[l];
    begin[l] = used; // new beginning of line
    debug_assert!(iend - ibeg <= fmem - used);
    for pos in ibeg..iend {
        index[used as usize] = index[pos as usize];
        value[used as usize] = value[pos as usize];
        used += 1;
    }
    end[l] = used;
    debug_assert!(fmem - used >= extra_space);
    used += extra_space;
    begin[n] = used; // beginning of unused space
    lu_list_move(line, 0, next, prev, nlines, None);
}

/// HiGHS: basiclu/lu_file.c::lu_file_compress (lines 113-144)
///
/// C `extra_space = stretch*(iend-ibeg) + pad` is a double truncated to
/// `lu_int`; `as LuInt` truncates toward zero the same way.
#[allow(clippy::too_many_arguments)] // C signature (BASICLU file API)
pub fn lu_file_compress(
    nlines: LuInt,
    begin: &mut [LuInt],
    end: &mut [LuInt],
    next: &[LuInt],
    index: &mut [LuInt],
    value: &mut [f64],
    stretch: f64,
    pad: LuInt,
) -> LuInt {
    let n = nlines as usize;
    let mut used: LuInt = 0;
    let mut extra_space: LuInt = 0;
    let mut nz: LuInt = 0;
    let mut i = next[n];
    while i < nlines {
        let iu = i as usize;
        let ibeg = begin[iu];
        let iend = end[iu];
        debug_assert!(ibeg >= used);
        used += extra_space;
        if used > ibeg {
            used = ibeg; // chop extra space added before
        }
        begin[iu] = used;
        for pos in ibeg..iend {
            index[used as usize] = index[pos as usize];
            value[used as usize] = value[pos as usize];
            used += 1;
        }
        end[iu] = used;
        extra_space = (stretch * ((iend - ibeg) as f64) + (pad as f64)) as LuInt;
        nz += iend - ibeg;
        i = next[iu];
    }
    debug_assert!(used <= begin[n]);
    used += extra_space;
    if used > begin[n] {
        used = begin[n]; // never use more space than before
    }
    begin[n] = used;
    nz
}

/// HiGHS: basiclu/lu_file.c::lu_file_diff (lines 164-184)
///
/// Only called inside `assert`s (`#ifndef NDEBUG`); use from `debug_assert!`.
pub fn lu_file_diff(
    nrow: LuInt,
    begin_row: &[LuInt],
    end_row: &[LuInt],
    begin_col: &[LuInt],
    end_col: &[LuInt],
    index: &[LuInt],
    value: Option<&[f64]>,
) -> LuInt {
    let mut ndiff: LuInt = 0;
    for i in 0..nrow {
        for pos in begin_row[i as usize]..end_row[i as usize] {
            let j = index[pos as usize] as usize;
            let end = end_col[j];
            let mut where_ = begin_col[j];
            while where_ < end && index[where_ as usize] != i {
                where_ += 1;
            }
            // C: `where == end || (value && value[pos] != value[where])`, short-circuit
            let missing = where_ == end
                || match value {
                    Some(v) => v[pos as usize] != v[where_ as usize],
                    None => false,
                };
            if missing {
                ndiff += 1;
            }
        }
    }
    ndiff
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-computed on a two-line file: line 0 holds 3 entries, line 1 holds 2,
    /// stretch 0.5 and pad 0. The extra space after line 0 is `0.5*3 + 0 = 1.5`,
    /// which truncates to 1 (C `lu_int` assignment), not 2.
    #[test]
    fn compress_truncates_extra_space_like_c() {
        let mut begin = vec![0, 10, 12];
        let mut end = vec![3, 12, 20];
        let next = vec![1, 2, 0];
        let mut index: Vec<LuInt> = vec![0; 20];
        let mut value: Vec<f64> = vec![0.0; 20];
        index[0..3].copy_from_slice(&[1, 2, 3]);
        value[0..3].copy_from_slice(&[1.0, 2.0, 3.0]);
        index[10..12].copy_from_slice(&[4, 5]);
        value[10..12].copy_from_slice(&[4.0, 5.0]);

        let nz = lu_file_compress(
            2, &mut begin, &mut end, &next, &mut index, &mut value, 0.5, 0,
        );

        assert_eq!(nz, 5);
        assert_eq!(begin, vec![0, 4, 7]);
        assert_eq!(end[..2], [3, 6]);
        assert_eq!(index[4..6], [4, 5]);
        assert_eq!(value[4..6], [4.0, 5.0]);
    }

    /// Empty file: heads of the memory list, unused space starts at 0.
    #[test]
    fn empty_sets_heads_and_links() {
        let mut begin = vec![9; 3];
        let mut end = vec![9; 3];
        let mut next = vec![9; 3];
        let mut prev = vec![9; 3];
        lu_file_empty(2, &mut begin, &mut end, &mut next, &mut prev, 10);

        assert_eq!(begin, vec![0, 0, 0]);
        assert_eq!(end, vec![0, 0, 10]);
        assert_eq!(next, vec![1, 2, 0]);
        assert_eq!(prev, vec![2, 0, 1]);
    }

    /// Line 0 (2 entries) is moved to the end of the used space and given room
    /// for 3 more entries; the unused-space pointer moves past that room.
    #[test]
    fn reappend_moves_line_to_end() {
        let mut begin = vec![0, 0, 2];
        let mut end = vec![2, 0, 10];
        let mut next = vec![2, 1, 0];
        let mut prev = vec![2, 1, 0];
        let mut index: Vec<LuInt> = vec![1, 2, 0, 0, 0, 0, 0, 0, 0, 0];
        let mut value: Vec<f64> = vec![1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        lu_file_reappend(
            0, 2, &mut begin, &mut end, &mut next, &mut prev, &mut index, &mut value, 3,
        );

        assert_eq!(begin, vec![2, 0, 7]);
        assert_eq!(end, vec![4, 0, 10]);
        assert_eq!(index[2..4], [1, 2]);
        assert_eq!(value[2..4], [1.0, 2.0]);
    }

    /// Row file of a 1x2 pattern, column file with one (0, val) entry per
    /// column. Values agree, so only a missing entry or a value mismatch counts.
    #[test]
    fn diff_counts_missing_and_mismatched_entries() {
        let begin_row = [0];
        let end_row = [2];
        let begin_col = [2, 3];
        let end_col = [3, 4];
        let mut index: Vec<LuInt> = vec![0, 1, 0, 0];
        let mut value: Vec<f64> = vec![1.0, 2.0, 1.0, 2.0];

        assert_eq!(
            lu_file_diff(
                1,
                &begin_row,
                &end_row,
                &begin_col,
                &end_col,
                &index,
                Some(&value)
            ),
            0
        );
        assert_eq!(
            lu_file_diff(1, &begin_row, &end_row, &begin_col, &end_col, &index, None),
            0
        );

        value[3] = 2.5;
        assert_eq!(
            lu_file_diff(
                1,
                &begin_row,
                &end_row,
                &begin_col,
                &end_col,
                &index,
                Some(&value)
            ),
            1
        );
        assert_eq!(
            lu_file_diff(1, &begin_row, &end_row, &begin_col, &end_col, &index, None),
            0
        );

        index[3] = 5;
        assert_eq!(
            lu_file_diff(1, &begin_row, &end_row, &begin_col, &end_col, &index, None),
            1
        );
    }
}
