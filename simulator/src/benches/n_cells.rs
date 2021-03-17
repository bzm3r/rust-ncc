use criterion::{criterion_group, criterion_main, Criterion};
use rust_ncc::{exp_setups, world};
use std::time::Duration;

fn n_cells() {
    let exp = exp_setups::n_cells::generate(Some(123829), 9);
    let mut w = world::World::new(exp, None, 10, 100);
    w.simulate(3.0 * 3600.0, false);
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("g0");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(120));
    group.bench_function("n_cells", |b| b.iter(|| n_cells()));
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
