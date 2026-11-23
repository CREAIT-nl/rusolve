//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! Ported from SCIP.

use super::*;
use crate::mip::symmetry_dynamic::lexred::lexicographic_reduction_add_permutation;
use crate::mip::symmetry_dynamic::orbitopal::{
    orbitopal_reduction_add_orbitope, orbitopal_reduction_default_column_ordering,
};
use std::collections::HashMap;

/// SCIP sym_rowcol.c: SCIP_ROWORDERING_BRANCHING
const S3_ROW_ORDERING_BRANCHING: HighsInt = 1;

/// SCIP 11.0.0 (32dd425): sym_rowcol.c::SCIPsymCompAddSymmetryHandling (lines 1040-1180)
/// and sym_lexorbred.c::addLexRed (lines 410-462)
///
/// Build S3SymData from HiGHS symmetries according to the D6 dispatch.
/// Returns None when there are no symmetries at all.
pub(crate) fn s3_build(
    symmetries: &HighsSymmetries,
    is_integral: &[bool],
) -> Option<Box<S3SymData>> {
    if symmetries.orbitopes.is_empty() && symmetries.num_perms == 0 {
        return None;
    }

    let mut data = S3SymData::default();
    let mut symvarmap: HashMap<HighsInt, HighsInt> = HashMap::new();
    let mut nsymvars: HighsInt = 0;

    // Dispatch over HiGHS orbitopes
    for (o, orbitope) in symmetries.orbitopes.iter().enumerate() {
        let nrows = orbitope.num_rows;
        let ncols = orbitope.row_length;

        // num_rows == 1: static (HiGHS keeps it for orbital_fixing)
        // SCIP sym_rowcol.c lines 1042-1071
        if nrows == 1 {
            data.static_orbitopes.push(o);
            continue;
        }

        // ncols == 2: LexRed permutation swapping column 0 and column 1 row by row
        // SCIP sym_rowcol.c lines 1084-1091
        if ncols == 2 {
            // permvars: 2*nrows matrix entries in column-major order
            let mut permvars: Vec<HighsInt> = Vec::with_capacity((2 * nrows) as usize);
            for j in 0..2 {
                for i in 0..nrows {
                    permvars.push(orbitope.matrix[(i + j * nrows) as usize]);
                }
            }

            // perm: perm[i + j*nrows] = i + (1-j)*nrows (swaps column 0 and column 1)
            let mut perm: Vec<HighsInt> = vec![0; (2 * nrows) as usize];
            for j in 0..2 {
                for i in 0..nrows {
                    perm[(i + j * nrows) as usize] = i + (1 - j) * nrows;
                }
            }

            let success = lexicographic_reduction_add_permutation(
                &mut data.lexdatas,
                &mut symvarmap,
                &mut nsymvars,
                &permvars,
                &perm,
                true,
            );
            if success {
                continue;
            }
            // fall through to dynamic if lexred failed
        }

        // num_set_packing_rows >= 3 and every row is set packing: static packing orbitope
        // SCIP sym_rowcol.c lines 1109-1146
        if orbitope.num_set_packing_rows >= 3 && orbitope.num_set_packing_rows == nrows {
            data.static_orbitopes.push(o);
            continue;
        }

        // Otherwise: dynamic orbitopal reduction
        // SCIP sym_rowcol.c lines 1160-1178
        let nelem = (nrows * ncols) as usize;
        let mut vars: Vec<HighsInt> = vec![0; nelem];
        let mut pos = 0usize;
        for i in 0..nrows {
            for j in 0..ncols {
                vars[pos] = orbitope.matrix[(i + j * nrows) as usize];
                pos += 1;
            }
        }

        let colordering = orbitopal_reduction_default_column_ordering();
        orbitopal_reduction_add_orbitope(
            &mut data,
            S3_ROW_ORDERING_BRANCHING,
            colordering,
            &vars,
            nrows,
            ncols,
            is_integral,
        );
    }

    // Components without an orbitope: add remaining permutations to LexRed
    // SCIP sym_lexorbred.c addLexRed lines 444-450 (`*success |= locsuccess` has no
    // reader here: the data is returned whenever symmetries exist)
    // PORT NOTE: HiGHS stores permutation images as column ids over
    // `permutation_columns`; `column_position` maps them to positions, which is the
    // form `lexicographic_reduction_add_permutation` takes (SCIP perm[] indices).
    let len = symmetries.permutation_columns.len();
    for p in 0..symmetries.num_perms {
        let base = p as usize * len;
        let mut perm: Vec<HighsInt> = Vec::with_capacity(len);
        for j in 0..len {
            let col_id = symmetries.permutations[base + j];
            let pos = symmetries.column_position[col_id as usize];
            debug_assert!(
                pos != -1,
                "permutation image {} not in permutation_columns",
                col_id
            );
            perm.push(pos);
        }
        lexicographic_reduction_add_permutation(
            &mut data.lexdatas,
            &mut symvarmap,
            &mut nsymvars,
            &symmetries.permutation_columns,
            &perm,
            true,
        );
    }

    Some(Box::new(data))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_orbitope(
        num_rows: HighsInt,
        row_length: HighsInt,
        num_set_packing_rows: HighsInt,
        matrix: Vec<HighsInt>,
    ) -> HighsOrbitopeMatrix {
        HighsOrbitopeMatrix {
            row_length,
            num_rows,
            num_set_packing_rows,
            row_is_set_packing: vec![HighsOrbitopeMatrix::kRowPacking; num_rows as usize],
            matrix,
            ..Default::default()
        }
    }

    #[test]
    fn s3_setup_none_without_symmetries() {
        let sym = HighsSymmetries::default();
        let is_integral = vec![true; 4];
        let result = s3_build(&sym, &is_integral);
        assert!(result.is_none());
    }

    #[test]
    fn s3_setup_single_row_orbitope_static() {
        let mut sym = HighsSymmetries::default();
        // 1 row, 3 cols: matrix = [0, 1, 2]
        sym.orbitopes.push(make_orbitope(1, 3, 1, vec![0, 1, 2]));
        sym.num_generators = 1;
        let is_integral = vec![true; 3];
        let result = s3_build(&sym, &is_integral).unwrap();
        assert_eq!(result.static_orbitopes.len(), 1);
        assert_eq!(result.static_orbitopes[0], 0);
        assert_eq!(result.orbitopes.len(), 0);
        assert_eq!(result.lexdatas.len(), 0);
    }

    #[test]
    fn s3_setup_two_column_orbitope_lexred() {
        let mut sym = HighsSymmetries::default();
        // 2 rows, 2 cols: matrix column-major [0,1, 2,3] => (0,0)=0,(1,0)=1,(0,1)=2,(1,1)=3
        sym.orbitopes.push(make_orbitope(2, 2, 2, vec![0, 1, 2, 3]));
        sym.num_generators = 1;
        let is_integral = vec![true; 4];
        let result = s3_build(&sym, &is_integral).unwrap();
        // ncols==2 should produce a lexdata
        assert_eq!(result.lexdatas.len(), 1);
        assert!(result.lexdatas[0].isdynamic);
        assert_eq!(result.orbitopes.len(), 0);
        assert_eq!(result.static_orbitopes.len(), 0);
    }

    #[test]
    fn s3_setup_full_orbitope_dynamic() {
        let mut sym = HighsSymmetries::default();
        // 3x3 orbitope, not all set packing (num_set_packing_rows=0)
        let matrix: Vec<HighsInt> = (0..9).collect();
        sym.orbitopes.push(make_orbitope(3, 3, 0, matrix));
        sym.num_generators = 1;
        let is_integral = vec![true; 9];
        let result = s3_build(&sym, &is_integral).unwrap();
        assert_eq!(result.orbitopes.len(), 1);
        assert_eq!(result.orbitopes[0].nrows, 3);
        assert_eq!(result.orbitopes[0].ncols, 3);
        // columnordering from orbitopal_reduction_default_column_ordering is 4 (MEDIAN)
        // rowordering is S3_ROW_ORDERING_BRANCHING = 1
        assert_eq!(result.orbitopes[0].rowordering, 1);
        assert_eq!(result.orbitopes[0].columnordering, 4);
        assert_eq!(result.static_orbitopes.len(), 0);
    }

    #[test]
    fn s3_setup_packing_orbitope_static() {
        let mut sym = HighsSymmetries::default();
        // 3 set-packing rows over 3 cols: num_set_packing_rows == nrows == 3
        let matrix: Vec<HighsInt> = (0..9).collect();
        sym.orbitopes.push(make_orbitope(3, 3, 3, matrix));
        sym.num_generators = 1;
        let is_integral = vec![true; 9];
        let result = s3_build(&sym, &is_integral).unwrap();
        assert_eq!(result.static_orbitopes.len(), 1);
        assert_eq!(result.static_orbitopes[0], 0);
        assert_eq!(result.orbitopes.len(), 0);
    }

    #[test]
    fn s3_setup_plain_permutation_lexred() {
        let mut sym = HighsSymmetries::default();
        // No orbitopes, one permutation acting on columns [0,1,2]
        // permutation: 0->1, 1->2, 2->0 (cyclic shift)
        sym.permutation_columns = vec![0, 1, 2];
        sym.column_position = vec![0, 1, 2];
        sym.permutations = vec![1, 2, 0];
        sym.num_perms = 1;
        sym.num_generators = 1;
        let is_integral = vec![true; 3];
        let result = s3_build(&sym, &is_integral).unwrap();
        assert_eq!(result.lexdatas.len(), 1);
        assert!(result.lexdatas[0].isdynamic);
        assert_eq!(result.orbitopes.len(), 0);
        assert_eq!(result.static_orbitopes.len(), 0);
    }
}
