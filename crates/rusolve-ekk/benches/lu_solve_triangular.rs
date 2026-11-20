//! BASICLU sparse triangular solve (symbolic DFS + numeric sweep), the hot loop of the IPX
//! crossover push phase. Synthetic lower-triangular factor with 3
//! entries below the diagonal per column, so a single-entry rhs reaches a few percent of m.

use criterion::{Criterion, criterion_group, criterion_main};
use rusolve_ekk::ipm::basiclu::lu_solve_symbolic::lu_solve_symbolic;
use rusolve_ekk::ipm::basiclu::lu_solve_triangular::lu_solve_triangular;

fn bench(c: &mut Criterion) {
    let m: usize = 20000;
    let mut seed: u64 = 12345;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) as usize
    };
    let mut begin = Vec::with_capacity(m + 1);
    let mut index: Vec<i32> = Vec::new();
    let mut value: Vec<f64> = Vec::new();
    for j in 0..m {
        begin.push(index.len() as i32);
        for _ in 0..3 {
            // rows close below the diagonal keep the reach large
            let i = j + 1 + next() % 40;
            if i < m {
                index.push(i as i32);
                value.push(0.1 + (next() % 100) as f64 * 0.001);
            }
        }
        index.push(-1); // sentinel-terminated column, as BASICLU's L
        value.push(0.0);
    }
    begin.push(index.len() as i32);

    let mut group = c.benchmark_group("lu_solve_triangular");
    group.sample_size(20);
    group.bench_function("m20000_reach", |b| {
        let mut lhs = vec![0.0f64; m];
        let mut pattern = vec![0i32; m];
        let mut pattern_symb = vec![0i32; m];
        let mut work1 = vec![0.0f64; m];
        let mut marked = vec![0i32; m];
        let mut marker = 0;
        b.iter(|| {
            marker += 1;
            let irhs = [0i32];
            let top = lu_solve_symbolic(
                m as i32,
                &begin,
                None,
                &index,
                1,
                &irhs,
                &mut pattern_symb,
                &mut work1,
                &mut marked,
                marker,
            );
            lhs[0] = 1.0;
            let mut flops = 0;
            let nz = lu_solve_triangular(
                m as i32 - top,
                &pattern_symb[top as usize..],
                &begin,
                None,
                &index,
                &value,
                None,
                0.0,
                &mut lhs,
                &mut pattern,
                &mut flops,
            );
            for k in 0..nz as usize {
                lhs[pattern[k] as usize] = 0.0;
            }
            std::hint::black_box(nz)
        })
    });
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
