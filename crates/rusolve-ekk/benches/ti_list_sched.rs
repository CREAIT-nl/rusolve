//! Time-indexed list scheduling (`mip_ti_list_sched`): one hill-climbing call on a synthetic
//! job shop, the hot loop of the heuristic.

use criterion::{Criterion, criterion_group, criterion_main};
use rusolve_ekk::mip::ti_precedence_separator::TiBenchInstance;

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("ti_list_sched");
    group.sample_size(10);
    for &(jobs, machines, horizon) in &[(6usize, 6usize, 36usize), (10, 5, 150)] {
        let mut inst = TiBenchInstance::new(jobs, machines, horizon);
        group.bench_function(format!("{jobs}x{machines}_h{horizon}"), |b| {
            b.iter(|| std::hint::black_box(inst.run()))
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
