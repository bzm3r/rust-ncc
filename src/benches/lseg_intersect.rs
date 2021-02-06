use criterion::{
    black_box, criterion_group, criterion_main, Criterion,
};
use rand::distributions::Uniform;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_ncc::math::geometry::{IntersectCalcResult, LineSeg2D};
use rust_ncc::math::v2d::V2D;
use std::time::Duration;

const SEED: u64 = 23891;

fn random_lseg(rng: &mut Pcg32, distrib: &Uniform<f64>) -> LineSeg2D {
    let a = V2D::new(rng.sample(&distrib), rng.sample(&distrib));
    let b = V2D::new(rng.sample(&distrib), rng.sample(&distrib));
    LineSeg2D::new(&a, &b)
}

fn random_lsegs(n: usize) -> Vec<LineSeg2D> {
    let mut rng = Pcg32::seed_from_u64(SEED);
    let uniform = Uniform::new(0.0_f64, 10.0);
    let lsegs =
        (0..n).map(|_| random_lseg(&mut rng, &uniform)).collect();
    lsegs
}

fn check_equivalence(lsegs: &[LineSeg2D]) {
    for (i, lseg_a) in lsegs.iter().enumerate() {
        'inner: for (_, lseg_b) in lsegs[(i + 1)..].iter().enumerate()
        {
            let calc_result = lseg_a.calc_lseg_intersect(&lseg_b);
            match (
                lseg_a.check_strict_lseg_intersect(&lseg_b),
                calc_result,
            ) {
                (true, IntersectCalcResult::Strict(_, _))
                | (true, IntersectCalcResult::Weak(_, _))
                | (false, _) => {
                    continue 'inner;
                }
                (true, _) => {
                    println!(
                        "lseg_a: {}, lseg_b: {}",
                        lseg_a, lseg_b
                    );
                    println!(
                        "check_lseg_intersect returns: {}",
                        lseg_a.check_strict_lseg_intersect(&lseg_b)
                    );
                    println!(
                        "calc_lseg_intersect returns: {}",
                        calc_result
                    );
                    panic!("unexpected difference between check_lseg_intersect and calc_lseg_intersect results!");
                }
            }
        }
    }
}

fn check_only(lsegs: &[LineSeg2D]) {
    for (i, lseg_a) in lsegs.iter().enumerate() {
        for lseg_b in lsegs[(i + 1)..].iter() {
            let _ = lseg_a.check_strict_lseg_intersect_v6(&lseg_b);
        }
    }
}

fn calc(lsegs: &[LineSeg2D]) {
    for (i, lseg_a) in lsegs.iter().enumerate() {
        for lseg_b in lsegs[(i + 1)..].iter() {
            let _ = lseg_a.calc_lseg_intersect(&lseg_b);
        }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let lsegs = random_lsegs(1000);
    let mut group = c.benchmark_group("g0");
    group.measurement_time(Duration::from_secs(40));
    group.bench_function("check_only", |b| {
        b.iter(|| check_only(black_box(&lsegs)))
    });
    group.bench_function("calc", |b| {
        b.iter(|| calc(black_box(&lsegs)))
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
