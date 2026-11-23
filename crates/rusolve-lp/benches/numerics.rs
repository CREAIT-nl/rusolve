//! Hot-path numerics benches (established with the first factorization slice).
//!
//! `cargo bench --package rusolve-lp --no-run` must compile; regression
//! numbers belong in `rusolve-bench`.

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use rusolve_core::CscMatrix;
use rusolve_lp::dense::DenseLu;
use rusolve_lp::rng::XorShift64;
use rusolve_lp::{DenseMatrix, PfiConfig, ProductFormFactorization};
use std::hint::black_box;

/// `y = A x` on a synthetic sparse matrix of ~10k nonzeros.
fn matvec_10k_nnz(c: &mut Criterion) {
    let n = 2000;
    let mut rng = XorShift64::new(20261021);
    let mut triplets = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let i = (rng.next_u64() % n as u64) as usize;
        let j = (rng.next_u64() % n as u64) as usize;
        triplets.push((i, j, rng.uniform(-1.0, 1.0)));
    }
    let a = CscMatrix::from_triplets(n, n, &triplets);
    let x: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut y = vec![0.0; n];

    c.bench_function("csc_matvec_10k_nnz", |b| b.iter(|| a.matvec(&x, &mut y)));
}

/// Triangular-ish `n x n` basis: dominant diagonal, dense lower band plus a
/// few random far entries — the shape simplex bases tend to have.
fn triangularish_basis(n: usize, seed: u64) -> DenseMatrix {
    let mut rng = XorShift64::new(seed);
    let mut a = DenseMatrix::zeros(n);
    for j in 0..n {
        for i in 0..n {
            let v = if i > j {
                rng.uniform(-0.5, 0.5)
            } else if i == j {
                rng.uniform(1.5, 2.5)
            } else {
                // Sparse upper part: ~2 far entries per column.
                0.0
            };
            a.set(i, j, v);
        }
        if j > 0 {
            for _ in 0..2 {
                let i = (rng.next_u64() % j as u64) as usize; // strictly above diagonal
                a.set(i, j, rng.uniform(-0.5, 0.5));
            }
        }
    }
    a
}

/// `solve_B` on a 500x500 triangular-ish basis: PFI (fresh, no etas) vs the
/// bare dense LU it wraps — isolates the product-form dispatch overhead.
fn solve_b_500(c: &mut Criterion) {
    let n = 500;
    let basis = triangularish_basis(n, 20261022);
    let lu = DenseLu::factor(&basis, 1e-12).expect("factorizable");
    let pfi = ProductFormFactorization::build(
        &basis,
        PfiConfig {
            pivot_tol: 1e-10,
            eta_limit: 1000,
            min_updates: 50,
            clock_scale: 1.0,
            engine: rusolve_lp::PfiEngine::Dense,
        },
    )
    .expect("buildable");
    let b: Vec<f64> = {
        let mut rng = XorShift64::new(7);
        (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect()
    };
    let mut x = vec![0.0; n];

    let mut g = c.benchmark_group("solve_b_500");
    g.bench_function("dense_lu", |ben| {
        ben.iter(|| lu.solve(black_box(&b), &mut x))
    });
    g.bench_function("pfi_no_etas", |ben| {
        ben.iter(|| pfi.solve_b(black_box(&b), &mut x))
    });
    g.finish();
}

/// Cost of one eta update as a function of current eta-file length (each
/// update embeds a full solve through the file, so this is the simplex
/// iteration-cost proxy that drives the refactorization threshold).
fn eta_update_cost(c: &mut Criterion) {
    let n = 500;
    let basis = triangularish_basis(n, 20261023);
    let cfg = PfiConfig {
        pivot_tol: 1e-10,
        eta_limit: usize::MAX,
        min_updates: 50,
        clock_scale: 1.0,
        engine: rusolve_lp::PfiEngine::Dense,
    };

    let mut g = c.benchmark_group("eta_update_cost");
    for &file_len in &[0usize, 10, 25, 50] {
        g.bench_with_input(format!("etas_{file_len:03}"), &file_len, |ben, &k| {
            ben.iter_batched(
                || {
                    let mut pfi = ProductFormFactorization::build(&basis, cfg).unwrap();
                    let mut rng = XorShift64::new(1000 + k as u64);
                    // Warm the file to exactly k etas with max-pivot updates.
                    for _ in 0..k {
                        let entering: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
                        let mut bbar = vec![0.0; n];
                        pfi.solve_b(&entering, &mut bbar);
                        let pos = (0..n)
                            .max_by(|i, j| bbar[*i].abs().total_cmp(&bbar[*j].abs()))
                            .unwrap();
                        pfi.update(pos, &entering, &bbar, 1.0).unwrap();
                    }
                    let next: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0, 1.0)).collect();
                    let mut bbar = vec![0.0; n];
                    pfi.solve_b(&next, &mut bbar);
                    let pos = (0..n)
                        .max_by(|i, j| bbar[*i].abs().total_cmp(&bbar[*j].abs()))
                        .unwrap();
                    (pfi, pos, next)
                },
                |(mut pfi, pos, next)| {
                    // `update` does not re-solve internally — the
                    // pivot's FTRAN feeds its b̄ from the caller. Measure
                    // solve + update together to keep this the per-pivot
                    // cost proxy it always was.
                    let mut bbar = vec![0.0; n];
                    pfi.solve_b(black_box(&next), &mut bbar);
                    pfi.update(black_box(pos), black_box(&next), black_box(&bbar), 1.0)
                        .unwrap();
                },
                BatchSize::LargeInput,
            );
        });
    }
    g.finish();
}

/// Factorization cost on a 500×500: dense workspace (DenseLu) vs the
/// sparse workspace (SparseLu). This guards the factorization cost that
/// made an early PFI engine swap time out on tuff. Shapes: random density-2% and the
/// triangular-ish basis-like matrix.
fn factor_500_sparse_vs_dense(c: &mut Criterion) {
    let n = 500;
    let mut rng = XorShift64::new(20261024);
    let random_cols: Vec<Vec<f64>> = (0..n)
        .map(|_| {
            (0..n)
                .map(|_| {
                    if rng.uniform(0.0, 1.0) < 0.02 {
                        rng.uniform(-10.0, 10.0)
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let tri = triangularish_basis(n, 20261025);
    let tri_cols: Vec<Vec<f64>> = (0..n).map(|j| tri.col(j).to_vec()).collect();
    // Basis-adjacent shape: SPARSE band + diagonal (bandwidth 8, ~3.4%
    // density, fill-controllable) — unlike the expander-random and the
    // structurally dense-lower triangularish shapes above, both of which
    // are adversarial for any sparse LU (no pivoting escapes an expander;
    // a dense lower triangle is dense work whatever the lists do).
    let mut rng2 = XorShift64::new(20261026);
    let band_cols: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            (0..n)
                .map(|i| {
                    let d = i as i64 - j as i64;
                    if d == 0 {
                        rng2.uniform(1.5, 2.5)
                    } else if d.abs() <= 8 && rng2.uniform(0.0, 1.0) < 0.5 {
                        rng2.uniform(-1.0, 1.0)
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();

    let mut g = c.benchmark_group("factor_500");
    g.bench_with_input("dense_random2pct", &random_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| DenseLu::factor(&DenseMatrix::from_columns(&cols), 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.bench_with_input("sparse_random2pct", &random_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| rusolve_lp::lu::SparseLu::factor(&cols, 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.bench_with_input("dense_triangularish", &tri_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| DenseLu::factor(&DenseMatrix::from_columns(&cols), 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.bench_with_input("sparse_triangularish", &tri_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| rusolve_lp::lu::SparseLu::factor(&cols, 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.bench_with_input("dense_band8", &band_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| DenseLu::factor(&DenseMatrix::from_columns(&cols), 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.bench_with_input("sparse_band8", &band_cols, |ben, cols| {
        ben.iter_batched(
            || cols.clone(),
            |cols| rusolve_lp::lu::SparseLu::factor(&cols, 1e-12).unwrap(),
            BatchSize::LargeInput,
        )
    });
    g.finish();
}

/// Materialization cost: the row-form front door
/// (`to_standard` — validation + the `from_lp` bounds hop + the CSC-only
/// construction) vs the caller-space CSC door (`to_standard_sparse` direct)
/// on the CFL-large shape (3,130 rows × 3,030 cols, nnz 12,030, 0.13 %
/// density — a facility-location shape whose per-node LP cost motivated the
/// sparse input path). Both build the same `Standard` (pinned by
/// `to_standard_delegates_to_the_csc_construction` in simplex.rs; the
/// load-bearing check is the Netlib A/B) — this isolates the extra
/// copy the row-form hop pays.
fn to_standard_materialize(c: &mut Criterion) {
    let (m, n, nnz) = (3130usize, 3030usize, 12_030usize);
    let mut rng = XorShift64::new(20261027);
    let mut triplets = Vec::with_capacity(nnz);
    for _ in 0..nnz {
        let i = (rng.next_u64() % m as u64) as usize;
        let j = (rng.next_u64() % n as u64) as usize;
        triplets.push((i, j, rng.uniform(0.25, 4.0)));
    }
    let a = CscMatrix::from_triplets(m, n, &triplets);
    let lp = rusolve_lp::Lp {
        num_col: n,
        num_row: m,
        a: a.clone(),
        rhs: vec![1.0; m],
        objective: vec![1.0; n],
        maximize: false,
        senses: vec!['G'; m],
        lower: vec![0.0; n],
        upper: vec![f64::INFINITY; n],
    };
    let sp = rusolve_core::SparseLp::from_lp(&lp);

    let mut g = c.benchmark_group("to_standard_materialize");
    // One construction path (the dense twin is gone); the bench
    // tracks the CSC build's cost.
    g.bench_function("sparse_cfl_shape", |ben| {
        ben.iter(|| {
            rusolve_lp::simplex::bench_hooks::standard_size(
                &rusolve_lp::simplex::bench_hooks::build_sparse(black_box(&sp)),
            )
        })
    });
    g.finish();
}

/// Large-m `factor_sparse` on a near-slack basis shaped like large
/// staircase LP bases (stormG2/pds: m ≈ 0.1–0.5M,
/// ~85 % unit columns, a small Markowitz kernel). `cold` = fresh workspace
/// plus flat transpose (the older build); `warm` = persistent [`FactorWs`]
/// fed the row mirror, with the previous factorization recycled — the PFI
/// refactor path. Both produce bit-identical factors (pinned by
/// `factor_ws_bit_identical` in lu.rs).
fn factor_sparse_large(c: &mut Criterion) {
    use rusolve_lp::lu::{FactorWs, SparseLu};
    let m = 200_000usize;
    let mut rng = XorShift64::new(20261028);
    // Structural "home" rows: every 7th row; the rest carry a slack.
    let is_home = |i: usize| i % 7 == 3;
    let homes: Vec<usize> = (0..m).filter(|&i| is_home(i)).collect();
    let mut cols: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
    for i in 0..m {
        if !is_home(i) {
            cols.push(vec![(i, 1.0)]);
            continue;
        }
        let mut c = vec![(i, rng.uniform(2.0, 3.0))];
        for _ in 0..3 {
            let r = (rng.next_u64() % m as u64) as usize;
            // Mostly slack rows (peeled); ~2 % land on another home row
            // (forms the kernel).
            let r = if rng.uniform(0.0, 1.0) < 0.02 {
                homes[(rng.next_u64() % homes.len() as u64) as usize]
            } else if is_home(r) {
                r + 1
            } else {
                r
            };
            if r < m && c.iter().all(|&(k, _)| k != r) {
                c.push((r, rng.uniform(-0.5, 0.5)));
            }
        }
        c.sort_by_key(|&(k, _)| k);
        cols.push(c);
    }
    let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); m];
    for (j, col) in cols.iter().enumerate() {
        for &(i, v) in col {
            rows[i].push((j as u32, v));
        }
    }
    let mut g = c.benchmark_group("factor_sparse_200k");
    g.sample_size(20);
    g.bench_function("cold_transpose", |b| {
        b.iter(|| black_box(SparseLu::factor_sparse(&cols, 0.1).unwrap()))
    });
    let mut ws = FactorWs::default();
    g.bench_function("warm_mirror_recycled", |b| {
        b.iter(|| {
            let lu = SparseLu::factor_sparse_ws(&cols, Some(&rows), 0.1, &mut ws).unwrap();
            black_box(lu.build_work());
            ws.recycle(Box::new(lu));
        })
    });
    g.finish();
}

criterion_group!(
    numerics,
    matvec_10k_nnz,
    solve_b_500,
    eta_update_cost,
    factor_500_sparse_vs_dense,
    factor_sparse_large,
    to_standard_materialize
);
criterion_main!(numerics);
