//! IPM hot path: numeric supernodal Cholesky of `A·Θ·Aᵀ + δI`
//! (symbolic analysis once, numeric factor per IPM iteration) plus the
//! triangular solves. Inputs: the node-arc incidence matrix of a 2-D grid
//! (normal matrix = grid Laplacian: ordering-sensitive fill) and, when
//! present, Netlib instances from benches/data/netlib after IPM formulation.
//! Also the dense kernels (packed GEMM, blocked panel Cholesky; throughput
//! in flops via `Throughput::Elements`), the fundamental-supernode partition
//! (`relax0`, no amalgamation) for A/B, and any MPS files listed in
//! `RUSOLVE_BENCH_MPS` (comma-separated paths), e.g. relaxed MIPLIB LPs.
//!
//! Check: `cargo bench --package rusolve-lp --bench ipm_chol --no-run`.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use rusolve_core::CscMatrix;
use rusolve_lp::ipm::chol::NormalChol;
use rusolve_lp::ipm::kernels::{self, GemmWork, Pivot};
use rusolve_lp::rng::XorShift64;
use std::hint::black_box;

fn grid_incidence(k: usize) -> CscMatrix {
    let id = |i: usize, j: usize| i * k + j;
    let mut trip = Vec::new();
    let mut col = 0;
    for i in 0..k {
        for j in 0..k {
            if j + 1 < k {
                trip.push((id(i, j), col, 1.0));
                trip.push((id(i, j + 1), col, -1.0));
                col += 1;
            }
            if i + 1 < k {
                trip.push((id(i, j), col, 1.0));
                trip.push((id(i + 1, j), col, -1.0));
                col += 1;
            }
        }
    }
    // one slack per node keeps M nonsingular
    for r in 0..k * k {
        trip.push((r, col, 1.0));
        col += 1;
    }
    CscMatrix::from_triplets(k * k, col, &trip)
}

fn bench_matrix(c: &mut Criterion, name: &str, a: &CscMatrix) {
    let mut rng = XorShift64::new(40);
    let theta: Vec<f64> = (0..a.ncol())
        .map(|_| 10f64.powf(rng.uniform(-6.0, 6.0)))
        .collect();
    let mut g = c.benchmark_group(format!("ipm_chol/{name}"));
    let mut ch0 = NormalChol::analyse_opts(a, usize::MAX, 0.0, 1).unwrap();
    g.throughput(Throughput::Elements(ch0.stats.flops as u64));
    g.bench_function("factor_relax0", |b| {
        b.iter(|| ch0.factor(black_box(&theta), 1e-10))
    });
    drop(ch0);
    let mut ch = NormalChol::analyse(a, usize::MAX).unwrap();
    eprintln!(
        "{name}: m {} nnz(A) {} nnz(M) {} nnz(L) {} supernodes {} (fundamental {}) flops {:.3e}",
        ch.stats.m,
        ch.stats.nnz_a,
        ch.stats.nnz_m,
        ch.stats.nnz_l,
        ch.stats.supernodes,
        ch.stats.supernodes_fund,
        ch.stats.flops
    );
    g.bench_function("factor", |b| b.iter(|| ch.factor(black_box(&theta), 1e-10)));
    g.throughput(Throughput::Elements(4 * ch.stats.nnz_l as u64));
    ch.factor(&theta, 1e-10);
    let rhs: Vec<f64> = (0..a.nrow()).map(|i| (i as f64).sin()).collect();
    let mut x = rhs.clone();
    let mut work = Vec::new();
    g.bench_function("solve", |b| {
        b.iter(|| {
            x.copy_from_slice(&rhs);
            ch.solve(black_box(&mut x), &mut work)
        })
    });
    g.finish();
}

fn bench_dense(c: &mut Criterion) {
    let mut rng = XorShift64::new(43);
    let mut g = c.benchmark_group("ipm_chol/dense");
    let (m, n, k) = (1024usize, 512usize, 256usize);
    let a: Vec<f64> = (0..m * k).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut cm = vec![0.0; m * n];
    let mut w = GemmWork::default();
    g.throughput(Throughput::Elements((2 * m * n * k) as u64));
    g.bench_function("gemm_nt_1024x512x256", |b| {
        b.iter(|| {
            kernels::gemm_nt(m, n, k, -1.0, &a, m, &a, m, &mut cm, m, false, &mut w);
        })
    });
    let nn = 1500usize;
    let bm: Vec<f64> = (0..nn * nn).map(|_| rng.uniform(-1.0, 1.0)).collect();
    let mut spd = vec![0.0; nn * nn];
    kernels::gemm_nt(
        nn, nn, nn, 1.0, &bm, nn, &bm, nn, &mut spd, nn, false, &mut w,
    );
    for i in 0..nn {
        spd[i + i * nn] += nn as f64;
    }
    let mut blk = spd.clone();
    let mut ws = vec![GemmWork::default()];
    g.throughput(Throughput::Elements((nn * nn * nn / 3) as u64));
    g.bench_function("panel_cholesky_1500", |b| {
        b.iter(|| {
            blk.copy_from_slice(&spd);
            kernels::panel_cholesky(
                &mut blk,
                nn,
                nn,
                64,
                1,
                &mut ws,
                |_, _| Pivot::Accept,
                || true,
            )
        })
    });
    g.finish();
}

fn ipm_chol(c: &mut Criterion) {
    bench_dense(c);
    bench_matrix(c, "grid60", &grid_incidence(60));
    if let Ok(list) = std::env::var("RUSOLVE_BENCH_MPS") {
        for path in list.split(',').filter(|p| !p.is_empty()) {
            let Ok(model) = rusolve_io::read_model_file(path) else {
                continue;
            };
            let Ok(sp) = model.to_sparse_lp() else {
                continue;
            };
            let Ok(lp) = rusolve_lp::ipm::model::formulate(&sp) else {
                continue;
            };
            let name = std::path::Path::new(path)
                .file_stem()
                .map_or("mps".to_string(), |s| s.to_string_lossy().into_owned());
            bench_matrix(c, &name, &lp.a);
        }
    }
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../../benches/data/netlib");
    for inst in ["scsd6", "stair", "scrs8"] {
        let path = format!("{dir}/{inst}.mps");
        let Ok(model) = rusolve_io::read_model_file(&path) else {
            continue;
        };
        let Ok(sp) = model.to_sparse_lp() else {
            continue;
        };
        let Ok(lp) = rusolve_lp::ipm::model::formulate(&sp) else {
            continue;
        };
        bench_matrix(c, inst, &lp.a);
    }
}

criterion_group!(benches, ipm_chol);
criterion_main!(benches);
