use rayon::prelude::*;
use std::time::Instant;
use dashmap::DashMap;

fn gen_xys_f32() -> Vec<(f32, f32)> {
    let nverts = 100*16;
    let mut xys: Vec<(f32, f32)> = Vec::with_capacity(nverts);

    for i in 0..nverts {
        let val = i as f32;
        xys.push((val, val));
    }

    xys
}

fn ix_2d(n: usize, num_pairs: usize, ix: usize) -> (usize, usize) {
    let q = num_pairs - ix - 1;
    let mut y = (0.5*(1.0 + (1.0 + 8.0*(q as f32)).sqrt()).floor()) as usize;
    if y > n - 1 {
        y -= 1;
    }
    let mut i = n - y - 1;
    let b = num_pairs - (n - i)*(n - i - 1)/2;
    let mut j = (ix - b) + (i + 1);
    if j > n - 1 {
        i += 1;
        j -= 1;
    }
    (i, j)

}

fn dist_calc_seq(xys: &[(f32, f32)]) -> Vec<f32> {
    let n = xys.len();
    let num_pairs = n * (n - 1) / 2;
    let mut dist_vec: Vec<f32> = Vec::with_capacity(n);

    (0..num_pairs).for_each(|ix| {
        let (i, j) = ix_2d(n, num_pairs, ix);
        let (x1, y1) = xys[i];
        let (x2, y2) = xys[j];
        dist_vec.push((x1 - x2).powi(2) + (y1 - y2).powi(2));
    });

    dist_vec
}

fn dist_calc_par(xys: &[(f32, f32)]) -> Vec<f32> {
    let nverts = xys.len();
    let n = nverts * (nverts - 1) / 2;
    let mut dist_vec: Vec<f32> = vec![f32::default(); n];
    dist_vec.par_iter_mut().enumerate().for_each(|(ix, e)| {
        let (i, j) = ix_2d(nverts, n, ix);
        let (x1, y1) = xys[i];
        let (x2, y2) = xys[j];
        *e = (x1 - x2).powi(2) + (y1 - y2).powi(2);
    });

    dist_vec
}

fn dist_calc_chmap(xys: &[(f32, f32)]) -> DashMap<(usize, usize), f32> {
    let nverts = xys.len();
    let n = nverts*(nverts - 1)/2;
    let cmap: DashMap<(usize, usize), f32> = DashMap::with_capacity(n);

    (0..nverts).into_par_iter().for_each(
        |i| {
            ((i + 1)..nverts).into_par_iter().for_each(
                |j| {
                    let (x1, y1) = xys[i];
                    let (x2, y2) = xys[j];
                    cmap.insert((i, j), (x1 - x2).powi(2) + (y1 - y2).powi(2));
                }
            )
        });

    cmap
}

fn main() {
    let xys = gen_xys_f32();
    let nverts = xys.len();
    let n = nverts*(nverts - 1)/2;

    let now = Instant::now();
    let rseq = dist_calc_seq(&xys);
    let tseq = now.elapsed().as_millis();
    assert_eq!(rseq.len(), n);

    let now = Instant::now();
    let rpar = dist_calc_par(&xys);
    let tpar = now.elapsed().as_millis();
    assert_eq!(rpar.len(), n);

    let now = Instant::now();
    let rchmap = dist_calc_chmap(&xys);
    let tchmap = now.elapsed().as_millis();

    println!(
        "results match (seq, par): {}",
            rseq.iter()
            .zip(rpar.iter())
            .all(|(a, b)| {
                (a - b).abs() < std::f32::EPSILON
            })
    );

    println!("results match (seq, chmap): {}", rseq.iter().enumerate().all(|(ix, a)| {
        let k = ix_2d(nverts, n, ix);
        let b = rchmap.get(&k).unwrap().clone();
        (a - b).abs() < std::f32::EPSILON
    }));

    println!("sequential: {} ms", tseq);
    println!("parallel: {} ms", tpar);
    println!("chmap: {} ms", tchmap);
}
//
// fn main() {
//     let nverts = 300*16;
//     let n = nverts*(nverts - 1)/2;
//     println!("num elements: {}", n);
//     for ix in 0..n {
//         let (i, j) = ix_2d(nverts, n, ix);
//         if i > nverts || j > nverts {
//             panic!("{} = ({}, {})", ix, i, j);
//         }
//     }
// }
