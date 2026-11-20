//! HiGHS source: ipm/ipx/guess_basis.h, guess_basis.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::control::{Control, textline};
use super::model::Model;
use super::sparse_matrix::SparseMatrix;
use super::sparse_utils::{augmenting_path, depth_first_search};
use super::utils::sortperm;
use super::*;
use crate::util::highs_utils::cpp_max;

/// HiGHS: ipx/guess_basis.cc::ComputePattern (static, lines 11-25)
fn compute_pattern(
    L: &SparseMatrix,
    AI: &SparseMatrix,
    j: Int,
    rownumber: &[Int],
    pattern: &mut [Int],
    pstack: &mut [Int],
    marked: &mut [Int],
    marker: Int,
) -> Int {
    let m = L.rows();
    let Lp = L.colptr();
    let Li = L.rowidx();
    let mut top = m;
    for p in AI.begin(j)..AI.end(j) {
        let i = AI.index(p);
        if marked[i as usize] != marker {
            top = depth_first_search(
                i,
                Lp,
                Li,
                Some(rownumber),
                top,
                pattern,
                marked,
                marker,
                pstack,
            );
        }
    }
    top
}

/// HiGHS: ipx/guess_basis.cc::ComputeValues (static, lines 30-60)
fn compute_values(
    L: &SparseMatrix,
    AI: &SparseMatrix,
    j: Int,
    rownumber: &[Int],
    pattern: &[Int],
    top: Int,
    lhs: &mut Vector,
) -> Int {
    let m = L.rows();
    let Lp = L.colptr();
    let Li = L.rowidx();
    let Lx = L.values();
    for t in top..m {
        lhs[pattern[t as usize] as usize] = 0.0;
    }
    for p in AI.begin(j)..AI.end(j) {
        // scatter RHS into lhs
        lhs[AI.index(p) as usize] = AI.value(p);
    }
    let mut lhsmax = 0.0; // maximum entry that has not been pivotal
    let mut imax: Int = -1; // corresponding position in lhs
    for t in top..m {
        let i = pattern[t as usize];
        let temp = lhs[i as usize];
        let k = rownumber[i as usize];
        if temp != 0.0 {
            if k >= 0 {
                // Row i has been pivotal; update lhs.
                for p in Lp[k as usize]..Lp[(k + 1) as usize] {
                    lhs[Li[p as usize] as usize] -= Lx[p as usize] * temp;
                }
            } else if temp.abs() > lhsmax {
                // Row i has not been pivotal; update lhsmax.
                lhsmax = temp.abs();
                imax = i;
            }
        }
    }
    imax
}

/// HiGHS: ipx/guess_basis.cc::ProcessFreeColumns (static, lines 66-107)
///
/// `int* active` is `&mut [i32]` (C `int`, not `Int`).
fn process_free_columns(
    control: &Control,
    model: &Model,
    weights: &[f64],
    basis: &mut Vec<Int>,
    rownumber: &mut [Int],
    active: &mut [i32],
    interrupt: &mut bool,
) {
    let m = model.rows();
    let n = model.cols();
    let AI = model.AI();

    // Workspace and output of left-looking incomplete LU factorization.
    // L is stored without the unit diagonal, row indices are (unpermuted) row
    // indices from AI. The upper triangular factor is not stored.
    let mut pattern: Vec<Int> = vec![0; m as usize];
    let mut work: Vec<Int> = vec![0; m as usize];
    let mut marked: Vec<Int> = vec![-1; m as usize];
    let mut lhs = Vector::new(m);
    let mut L = SparseMatrix::with_dims(m, 0);

    let mut num_free: Int = 0;
    for j in 0..(n + m) {
        if weights[j as usize] != f64::INFINITY {
            continue;
        }
        let top = compute_pattern(
            &L,
            AI,
            j,
            rownumber,
            &mut pattern,
            &mut work,
            &mut marked,
            j,
        );
        let imax = compute_values(&L, AI, j, rownumber, &pattern, top, &mut lhs);
        let pivot = if imax >= 0 { lhs[imax as usize] } else { 0.0 };
        if pivot.abs() > kLuDependencyTol {
            debug_assert!(rownumber[imax as usize] == -1);
            rownumber[imax as usize] = basis.len() as Int;
            basis.push(j);
            // Add new column to L. Discard off-diagonal entries that were not
            // in the pattern of AI[:,j].
            for t in AI.begin(j)..AI.end(j) {
                let i = AI.index(t);
                if rownumber[i as usize] < 0 && lhs[i as usize] != 0.0 {
                    L.push_back(i, lhs[i as usize] / pivot);
                }
            }
            L.add_column();
            num_free += 1;
        }
        active[j as usize] = 0;
        // rusolve BASICLU time limit (not in the HiGHS reference): interrupt check at the end of each pass
        *interrupt = control.interrupt_check(-1) != 0;
        if *interrupt {
            return;
        }
    }
    control.debug(
        1,
        &format!(
            "{}{}\n",
            textline("Number of free variables in starting basis:"),
            num_free
        ),
    );
}

/// HiGHS: ipx/guess_basis.cc::ProcessSingletons (static, lines 111-145)
///
/// Adds singleton columns to the basis if their entry is maximum up to a factor
/// 2.0 in its row.
fn process_singletons(
    control: &Control,
    model: &Model,
    weights: &[f64],
    basis: &mut Vec<Int>,
    rownumber: &mut [Int],
    active: &mut [i32],
    interrupt: &mut bool,
) {
    let m = model.rows();
    let AI = model.AI();
    let AT = model.AIt();
    let mut num_singletons: Int = 0;
    for i in 0..m {
        if rownumber[i as usize] >= 0 {
            continue;
        }
        let mut rowmax = 0.0;
        let mut max_singleton = 0.0;
        let mut jsingleton: Int = -1;
        for p in AT.begin(i)..AT.end(i) {
            let j = AT.index(p);
            if active[j as usize] == 0 {
                continue;
            }
            let a = AT.value(p).abs() * weights[j as usize];
            rowmax = cpp_max(rowmax, a);
            if a > max_singleton && AI.end(j) == AI.begin(j) + 1 {
                max_singleton = a;
                jsingleton = j;
            }
        }
        if max_singleton > 0.0 && max_singleton >= 0.5 * rowmax {
            rownumber[i as usize] = basis.len() as Int;
            basis.push(jsingleton);
            active[jsingleton as usize] = 0;
            num_singletons += 1;
        }
        // rusolve BASICLU time limit (not in the HiGHS reference): interrupt check at the end of each pass
        *interrupt = control.interrupt_check(-1) != 0;
        if *interrupt {
            return;
        }
    }
    control.debug(
        1,
        &format!(
            "{}{}\n",
            textline("Number of singletons in starting basis:"),
            num_singletons
        ),
    );
}

/// HiGHS: ipx/guess_basis.cc::ProcessRemaining (static, lines 150-203)
///
/// Computes a matching of active columns with rows that have not been pivotal.
/// Each matched column is added to the basis. Columns are processed in
/// decreasing order of their weight.
fn process_remaining(
    control: &Control,
    model: &Model,
    weights: &[f64],
    basis: &mut Vec<Int>,
    rownumber: &mut [Int],
    active: &mut [i32],
    interrupt: &mut bool,
) {
    let m = model.rows();
    let n = model.cols();
    let AI = model.AI();
    let colperm = sortperm(n + m, Some(weights), true);

    // jmatch[i] == -1  if row i is unmatched but eligible for matching
    // jmatch[i] == -2  if row i is excluded from being matched
    // jmatch[i] == j >= 0 if row i is matched to column j
    let mut jmatch: Vec<Int> = vec![-1; m as usize];
    for i in 0..m {
        if rownumber[i as usize] >= 0 {
            jmatch[i as usize] = -2;
        }
    }
    let mut marked: Vec<Int> = vec![-1; (n + m) as usize];
    let Ap = AI.colptr();
    let Ai = AI.rowidx();
    let mut cheap: Vec<Int> = Ap[..(n + m) as usize].to_vec();
    let mut work: Vec<Int> = vec![0; m as usize];
    let mut work2: Vec<Int> = vec![0; (m + 1) as usize];
    let mut work3: Vec<Int> = vec![0; (m + 1) as usize];
    let mut num_matched: Int = 0;
    let mut num_failed: Int = 0;

    for &j in &colperm {
        if active[j as usize] == 0 {
            continue;
        }
        if weights[j as usize] == 0.0 {
            break;
        }
        let matched = augmenting_path(
            j,
            Ap,
            Ai,
            &mut jmatch,
            &mut cheap,
            &mut marked,
            &mut work,
            &mut work2,
            &mut work3,
        );
        if matched {
            basis.push(j);
            num_matched += 1;
        } else {
            num_failed += 1;
        }
        if num_failed >= 10 * (m - basis.len() as Int) {
            break;
        }
        // rusolve BASICLU time limit (not in the HiGHS reference): interrupt check at the end of each pass
        *interrupt = control.interrupt_check(-1) != 0;
        if *interrupt {
            return;
        }
    }
    for i in 0..m {
        if jmatch[i as usize] >= 0 {
            // rownumber[i] = m indicates that row i was matched
            debug_assert!(rownumber[i as usize] < 0);
            rownumber[i as usize] = m;
        }
    }
    control.debug(
        1,
        &format!(
            "{}{}\n{}{}\n",
            textline("Number of other columns matched:"),
            num_matched,
            textline("Number of other columns failed:"),
            num_failed
        ),
    );
}

/// HiGHS: ipx/guess_basis.cc::GuessBasis (lines 205-231)
pub fn guess_basis(
    control: &Control,
    model: &Model,
    colweights: &[f64],
    interrupt: &mut bool,
) -> Vec<Int> {
    let m = model.rows();
    let n = model.cols();

    // basis starts empty and is filled one index at a time. rownumber[i] >= 0
    // iff row i was pivot row when a column was added to the basis. (The
    // specific value has a different meaning in each method.) A column is
    // "active" if it is eligible for being added to the basis.
    let mut basis: Vec<Int> = Vec::new();
    let mut rownumber: Vec<Int> = vec![-1; m as usize];
    let mut active: Vec<i32> = vec![1; (n + m) as usize];

    process_free_columns(
        control,
        model,
        colweights,
        &mut basis,
        &mut rownumber,
        &mut active,
        interrupt,
    );
    if *interrupt {
        return basis;
    }
    process_singletons(
        control,
        model,
        colweights,
        &mut basis,
        &mut rownumber,
        &mut active,
        interrupt,
    );
    if *interrupt {
        return basis;
    }
    process_remaining(
        control,
        model,
        colweights,
        &mut basis,
        &mut rownumber,
        &mut active,
        interrupt,
    );
    if *interrupt {
        return basis;
    }

    // Complete basis with unit columns.
    for i in 0..m {
        if rownumber[i as usize] < 0 {
            basis.push(n + i);
        }
    }
    debug_assert!(basis.len() as Int == m);
    basis
}
