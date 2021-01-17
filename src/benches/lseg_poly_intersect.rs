use criterion::{
    black_box, criterion_group, criterion_main, Criterion,
};
use rand::distributions::Uniform;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;
use rust_ncc::math::geometry::{LineSeg2D, Poly};
use rust_ncc::math::v2d::V2D;
use rust_ncc::NVERTS;
use std::f32::consts::PI;
use std::time::Duration;

const SEED: u64 = 23891;

fn star_poly(rng: &mut Pcg32, translation: V2D) -> Poly {
    let mut points = [V2D::zeros(); NVERTS];
    let rs: Vec<f32> = rng
        .sample_iter(&Uniform::from(0.25_f32..1.0))
        .take(NVERTS)
        .collect();
    let thetas: Vec<f32> = (0..NVERTS)
        .map(|k| (k as f32 / NVERTS as f32) * 2.0 * PI)
        .collect();
    points.iter_mut().enumerate().for_each(|(i, p)| {
        *p = V2D::new(
            rs[i] * f32::cos(thetas[i]),
            rs[i] * f32::sin(thetas[i]),
        ) + translation
    });
    Poly::from_points(&points)
}

fn star_polys(n: usize) -> Vec<Poly> {
    let mut rng = Pcg32::seed_from_u64(SEED);
    let mut polys: Vec<Poly> = Vec::with_capacity(n);
    let uniform = Uniform::from(0.0_f32..10.0);
    while polys.len() < n {
        let translation: Vec<f32> =
            (&mut rng).sample_iter(&uniform).take(2).collect();
        let translation = V2D::new(translation[0], translation[1]);
        let new_poly = star_poly(&mut rng, translation);
        if !polys
            .iter()
            .any(|poly| new_poly.bbox.intersects(&poly.bbox))
        {
            polys.push(new_poly);
        } else {
            break;
        }
    }
    polys
}

fn gen_test_lsegs(polys: &[Poly]) -> Vec<LineSeg2D> {
    let num_polys = polys.len();
    let mut lsegs = Vec::with_capacity(
        (polys.len() * polys.len() - 1) * (NVERTS as usize).pow(2)
            / 2,
    );
    for m in 0..num_polys {
        let poly_a = &polys[m];
        for n in (m + 1)..num_polys {
            let poly_b = &polys[n];
            for i in 0..NVERTS {
                for j in 0..NVERTS {
                    lsegs.push(LineSeg2D::new(
                        &poly_a.verts[i],
                        &poly_b.verts[j],
                    ));
                }
            }
        }
    }
    lsegs
}

fn lseg_intersects_poly_check(lsegs: &[LineSeg2D], polys: &[Poly]) {
    for lseg in lsegs {
        for poly in polys {
            let _ = lseg.check_poly_intersect(poly);
        }
    }
}

fn old_lseg_intersects_poly_check(
    lsegs: &[LineSeg2D],
    polys: &[Poly],
) {
    for lseg in lsegs {
        for poly in polys {
            let _ = lseg.old_check_poly_intersect(poly);
        }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let polys = star_polys(1000);
    let lsegs = gen_test_lsegs(&polys);
    let mut group = c.benchmark_group("g0");
    group.measurement_time(Duration::from_secs(15));
    group.bench_function("lseg_intersects_poly_check", |b| {
        b.iter(|| {
            lseg_intersects_poly_check(
                black_box(&lsegs),
                black_box(&polys),
            )
        })
    });
    group.bench_function("old_lseg_intersects_poly_check", |b| {
        b.iter(|| {
            old_lseg_intersects_poly_check(
                black_box(&lsegs),
                black_box(&polys),
            )
        })
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
