//! HiGHS source: ipm/ipx/symbolic_invert.h, symbolic_invert.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! TRAP: `RandomPermute` uses libc++ `std::default_random_engine`
//! (`minstd_rand`, seed 1) and libc++ `uniform_int_distribution<Int>`;
//! both are reproduced below. It runs on
//! every `SymbolicInvert` (via `Basis::SymbolicInvert`, `DensityInverse`,
//! ipm.cc:886 `Debug(4)`, and `LpSolver::SymbolicInvert`), not on the
//! default solve path. Checked bit for bit against the oracle
//! (`basis_ops.random_permute` in tests/fixtures/ipx/linalg).

use super::model::Model;
use super::sparse_matrix::{SparseMatrix, copy_columns, transpose};
use super::sparse_utils::{augmenting_path, depth_first_search};
use super::*;

/// libc++ `std::minstd_rand` (`linear_congruential_engine<uint_fast32_t,
/// 48271, 0, 2147483647>`, default seed 1): `min() == 1`, `max() ==
/// 2147483646`. `uint_fast32_t` is `uint32_t` on Apple arm64.
struct MinstdRand {
    x: u32,
}

impl MinstdRand {
    const MIN: u32 = 1;
    const MAX: u32 = 2_147_483_646;

    fn new() -> Self {
        MinstdRand { x: 1 }
    }

    fn next(&mut self) -> u32 {
        self.x = ((u64::from(self.x) * 48271) % 2_147_483_647) as u32;
        self.x
    }
}

/// libc++ `uniform_int_distribution<int>::operator()(g, {a, b})` with the
/// `__independent_bits_engine<minstd_rand, uint32_t>` it builds on every
/// call (libc++ `<__random/uniform_int_distribution.h>`). Only the branch
/// reachable for `b - a + 1 <= 2^30` is taken here (IPX dimensions are
/// `Int`, far below), asserted below.
fn uniform_int(g: &mut MinstdRand, a: Int, b: Int) -> Int {
    const DT: u32 = u32::BITS; // numeric_limits<uint32_t>::digits
    let rp = (b as u32).wrapping_sub(a as u32).wrapping_add(1);
    if rp == 1 {
        return a;
    }
    let mut w = DT - rp.leading_zeros() - 1;
    if rp & (u32::MAX >> (DT - w)) != 0 {
        w += 1;
    }
    // __independent_bits_engine(g, w): engine range R = max - min + 1 =
    // 2^31 - 2, m = floor(log2 R) = 30, so n = 1 and w0 = w for w <= 30.
    debug_assert!((1..=30).contains(&w));
    let r: u32 = MinstdRand::MAX - MinstdRand::MIN + 1;
    let y0: u32 = (r >> w) << w;
    let mask0: u32 = u32::MAX >> (DT - w);
    loop {
        let mut u;
        loop {
            u = g.next() - MinstdRand::MIN;
            if u < y0 {
                break;
            }
        }
        let s = u & mask0;
        if s < rp {
            return (s.wrapping_add(a as u32)) as Int;
        }
    }
}

/// HiGHS: ipx/symbolic_invert.cc::RandomPermute (static, lines 10-18)
fn random_permute(basis: &[Int]) -> Vec<Int> {
    let m = basis.len() as Int;
    let mut re = MinstdRand::new();
    let mut permuted_basis: Vec<Int> = basis.to_vec();
    for k in 0..m {
        let d = uniform_int(&mut re, 0, m - 1);
        permuted_basis.swap(k as usize, d as usize);
    }
    permuted_basis
}

/// Test hook for the golden replay of the file-static `RandomPermute`
/// (the oracle compiles the same code against the same libc++).
#[doc(hidden)]
pub fn random_permute_for_oracle(basis: &[Int]) -> Vec<Int> {
    random_permute(basis)
}

/// HiGHS: ipx/symbolic_invert.cc::Matching (static, lines 20-52)
fn matching(model: &Model, basis: &[Int]) -> Vec<Int> {
    let m = model.rows();
    let n = model.cols();
    let AI = model.AI();
    let Ap = AI.colptr();
    let Ai = AI.rowidx();
    let mut jmatch: Vec<Int> = vec![-1; m as usize];
    let mut cheap: Vec<Int> = Ap[..(n + m) as usize].to_vec();
    let mut marked: Vec<Int> = vec![-1; (n + m) as usize];
    let mut work: Vec<Int> = vec![0; m as usize];
    let mut work2: Vec<Int> = vec![0; (m + 1) as usize];
    let mut work3: Vec<Int> = vec![0; (m + 1) as usize];

    // Match singletons first, since then they are all cheap matches.
    for &j in basis {
        if Ap[(j + 1) as usize] == Ap[j as usize] + 1 {
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
            debug_assert!(matched);
        }
    }
    for &j in basis {
        if Ap[(j + 1) as usize] != Ap[j as usize] + 1 {
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
            debug_assert!(matched);
        }
    }
    jmatch
}

/// HiGHS: ipx/symbolic_invert.cc::Blockperm (static, lines 55-92)
fn blockperm(AI: &SparseMatrix, jmatch: &[Int], BT: &SparseMatrix) -> Vec<Vec<Int>> {
    let m = AI.rows();
    let Ap = AI.colptr();
    let Ai = AI.rowidx();
    let BTp = BT.colptr();
    let BTi = BT.rowidx();

    // Put row indices into istack, working from the bottom of the stack (m-1)
    // to the top (0), so that indices whose depth-first search finished first
    // appear further down the stack.
    let mut istack: Vec<Int> = vec![0; m as usize];
    let mut marked: Vec<Int> = vec![0; m as usize];
    let mut work: Vec<Int> = vec![0; m as usize];
    let mut top = m;
    for i in 0..m {
        if marked[i as usize] != 1 {
            top = depth_first_search(
                i,
                Ap,
                Ai,
                Some(jmatch),
                top,
                &mut istack,
                &mut marked,
                1,
                &mut work,
            );
        }
    }
    debug_assert!(top == 0);

    let mut blocks: Vec<Vec<Int>> = Vec::new();
    let mut rowperm: Vec<Int> = vec![0; m as usize];
    top = m;
    for &i in istack.iter() {
        if marked[i as usize] != 2 {
            // Node i is not part of any block yet.
            let end = top;
            top = depth_first_search(
                i,
                BTp,
                BTi,
                None,
                top,
                &mut rowperm,
                &mut marked,
                2,
                &mut work,
            );
            debug_assert!(top < end);
            blocks.push(rowperm[top as usize..end as usize].to_vec());
        }
    }
    debug_assert!(top == 0);
    blocks.reverse();
    blocks
}

/// HiGHS: ipx/symbolic_invert.cc::CoarsenedGraph (static, lines 95-126)
fn coarsened_graph(BT: &SparseMatrix, blocks: &[Vec<Int>]) -> SparseMatrix {
    let m = BT.rows();
    let nb = blocks.len() as Int;

    let mut map2block: Vec<Int> = vec![-1; m as usize];
    for b in 0..nb {
        for &i in &blocks[b as usize] {
            debug_assert!(map2block[i as usize] == -1);
            map2block[i as usize] = b;
        }
    }
    for i in 0..m {
        debug_assert!(map2block[i as usize] >= 0);
    }

    let mut cmat = SparseMatrix::with_dims(nb, 0);
    let mut marked: Vec<Int> = vec![-1; m as usize];
    for k in 0..nb {
        for &i in &blocks[k as usize] {
            for p in BT.begin(i)..BT.end(i) {
                let b = map2block[BT.index(p) as usize];
                if marked[b as usize] != k {
                    marked[b as usize] = k;
                    cmat.push_back(b, 1.0);
                    debug_assert!(b >= k);
                }
            }
        }
        cmat.add_column();
    }
    cmat
}

/// HiGHS: ipx/symbolic_invert.cc::SymbolicInvert (lines 128-181)
///
/// `rowcounts` / `colcounts` may be `nullptr`.
pub fn symbolic_invert(
    model: &Model,
    basis: &[Int],
    rowcounts: Option<&mut [Int]>,
    colcounts: Option<&mut [Int]>,
) {
    let AI = model.AI();
    let m = AI.rows();
    debug_assert!(basis.len() as Int == m);

    // jmatch is a permutation of basis such that B = AI[:,jmatch] has a
    // zero-free diagonal.
    let jmatch = matching(model, &random_permute(basis));
    let BT = transpose(&copy_columns(AI, &jmatch));

    let blocks = blockperm(AI, &jmatch, &BT);
    let mut cmat = coarsened_graph(&BT, &blocks);

    let nb = blocks.len() as Int;
    let mut stack: Vec<Int> = vec![0; nb as usize];
    let mut marked: Vec<Int> = vec![0; nb as usize];
    let mut work: Vec<Int> = vec![0; nb as usize];

    if let Some(rowcounts) = rowcounts {
        let mut jcount: Vec<Int> = vec![-1; AI.cols() as usize];
        marked.fill(-1);
        let Cp = cmat.colptr();
        let Ci = cmat.rowidx();
        for b in 0..nb {
            let top =
                depth_first_search(b, Cp, Ci, None, nb, &mut stack, &mut marked, b, &mut work);
            let mut nz: Int = 0;
            for t in top..nb {
                nz += blocks[stack[t as usize] as usize].len() as Int;
            }
            for &i in &blocks[b as usize] {
                jcount[jmatch[i as usize] as usize] = nz;
            }
        }
        for p in 0..m {
            rowcounts[p as usize] = jcount[basis[p as usize] as usize];
            debug_assert!(rowcounts[p as usize] >= 0);
        }
    }

    if let Some(colcounts) = colcounts {
        cmat = transpose(&cmat);
        marked.fill(-1);
        let Cp = cmat.colptr();
        let Ci = cmat.rowidx();
        for b in 0..nb {
            let top =
                depth_first_search(b, Cp, Ci, None, nb, &mut stack, &mut marked, b, &mut work);
            let mut nz: Int = 0;
            for t in top..nb {
                nz += blocks[stack[t as usize] as usize].len() as Int;
            }
            for &i in &blocks[b as usize] {
                colcounts[i as usize] = nz;
            }
        }
    }
}
