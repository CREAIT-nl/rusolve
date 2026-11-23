//! Builder cost per nonzero (docs/design/python-api.md P3). Throughput is reported in nnz, so
//! criterion's time / element is ns per nonzero.
use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use rusolve_model::{Fill, FillBool, ModelBuilder};
use std::hint::black_box;

const NCOLS: usize = 100_000;

fn model(ncols: usize) -> ModelBuilder {
    let mut b = ModelBuilder::new();
    b.add_cols(
        ncols,
        Fill::Scalar(1.0),
        Fill::Scalar(0.0),
        Fill::Scalar(f64::INFINITY),
        FillBool::Scalar(false),
    )
    .unwrap();
    b
}

/// Deterministic CSR with `nnz_per_row` distinct columns per row.
fn csr(nrows: usize, per_row: usize) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    let mut indptr = Vec::with_capacity(nrows + 1);
    let mut indices = Vec::with_capacity(nrows * per_row);
    let mut data = Vec::with_capacity(nrows * per_row);
    let mut x: u64 = 0x9E3779B97F4A7C15;
    indptr.push(0);
    for _ in 0..nrows {
        let start = indices.len();
        while indices.len() - start < per_row {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let j = (x % NCOLS as u64) as i32;
            if !indices[start..].contains(&j) {
                indices.push(j);
                data.push(1.0 + (x % 7) as f64);
            }
        }
        indptr.push(indices.len() as i32);
    }
    (indptr, indices, data)
}

fn bench(c: &mut Criterion) {
    let mut g = c.benchmark_group("builder");
    for per_row in [2usize, 50] {
        let nrows = 1_000_000 / per_row / 4; // 250k nnz per iteration
        let (indptr, indices, data) = csr(nrows, per_row);
        g.throughput(Throughput::Elements(indices.len() as u64));
        g.bench_function(format!("push_finish_{per_row}nnz"), |bch| {
            bch.iter_batched(
                || model(NCOLS),
                |mut b| {
                    for r in 0..nrows {
                        for p in indptr[r] as usize..indptr[r + 1] as usize {
                            b.push(indices[p] as i64, data[p]).unwrap();
                        }
                        b.finish_row(-1.0, 1.0).unwrap();
                    }
                    black_box(b)
                },
                BatchSize::LargeInput,
            )
        });
    }
    let (indptr, indices, data) = csr(1_000_000 / 10, 10);
    g.throughput(Throughput::Elements(indices.len() as u64));
    g.sample_size(20);
    g.bench_function("add_rows_csr_1M", |bch| {
        bch.iter_batched(
            || model(NCOLS),
            |mut b| {
                b.add_rows_csr(
                    &indptr,
                    &indices,
                    &data,
                    None,
                    Fill::Scalar(-1.0),
                    Fill::Scalar(1.0),
                )
                .unwrap();
                black_box(b)
            },
            BatchSize::LargeInput,
        )
    });
    g.bench_function("take_into_lp_1M", |bch| {
        bch.iter_batched(
            || {
                let mut b = model(NCOLS);
                b.add_rows_csr(
                    &indptr,
                    &indices,
                    &data,
                    None,
                    Fill::Scalar(-1.0),
                    Fill::Scalar(1.0),
                )
                .unwrap();
                b
            },
            |mut b| black_box(b.take().into_lp()),
            BatchSize::LargeInput,
        )
    });
    g.bench_function("csr_take_into_lp_1M_total", |bch| {
        bch.iter_batched(
            || model(NCOLS),
            |mut b| {
                b.add_rows_csr(
                    &indptr,
                    &indices,
                    &data,
                    None,
                    Fill::Scalar(-1.0),
                    Fill::Scalar(1.0),
                )
                .unwrap();
                black_box(b.take().into_lp())
            },
            BatchSize::LargeInput,
        )
    });
    g.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
