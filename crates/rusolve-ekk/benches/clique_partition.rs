//! Clique partition of the objective binaries (`setup_clique_partition`), the hot spot of the
//! FERC unit-commitment setup: neighbourhood scan against the cliques of the seed
//! (`h16_hash_common` on) and one hash lookup per candidate (off).

use criterion::{Criterion, criterion_group, criterion_main};
use rusolve_ekk::mip::clique_table::{HighsCliqueTable, HighsCliqueTable_CliqueVar};

/// `n` binaries, size-two cliques (i, i+1) and (i, i+7): the sparse adjacency of a unit
/// commitment start/stop structure.
fn table(n: i32, exact: bool) -> HighsCliqueTable {
    let mut t = HighsCliqueTable::new(n);
    t.h16_hash_common = exact;
    for i in 0..n {
        for d in [1, 7] {
            if i + d < n {
                let c = [
                    HighsCliqueTable_CliqueVar::new(i, 0),
                    HighsCliqueTable_CliqueVar::new(i + d, 0),
                ];
                t.do_add_clique(&c, false, i32::MAX);
            }
        }
    }
    t
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("clique_partition");
    group.sample_size(10);
    let n = 6000;
    let cost: Vec<f64> = (0..n).map(|i| 1.0 + (i % 13) as f64).collect();
    for (name, exact) in [("scan", false), ("cliques", true)] {
        let base = table(n, exact);
        group.bench_function(name, |b| {
            b.iter(|| {
                let mut t = base.clone();
                let mut vars: Vec<_> = (0..n)
                    .map(|i| HighsCliqueTable_CliqueVar::new(i, 0))
                    .collect();
                let mut start = Vec::new();
                t.clique_partition_2(&cost, &mut vars, &mut start);
                std::hint::black_box(start.len())
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
