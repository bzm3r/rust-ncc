use criterion::{criterion_group, criterion_main, Criterion};
use rust_ncc::{experiments, world};
use std::time::Duration;

fn n_cells() {
    let exp = experiments::n_cells::generate(Some(123829), 9);
    let mut w = world::World::new(exp, None, 10);
    w.simulate(3.0 * 3600.0);
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("g0");
    group.measurement_time(Duration::from_secs(120));
    group.bench_function("n_cells", |b| b.iter(|| n_cells()));
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
