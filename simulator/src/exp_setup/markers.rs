use crate::math::radians::{Radians, RAD_PI};
use crate::NVERTS;

pub fn mark_between_angles(bounds: (Radians, Radians)) -> [bool; NVERTS] {
    let mut r = [false; NVERTS];
    let (b0, b1) = bounds;
    (0..NVERTS).for_each(|vi| {
        let va = (2.0 * vi as f64 / NVERTS as f64) * (*RAD_PI);
        if va.between(b0, b1) {
            r[vi] = true;
        }
    });
    r
}

pub const ALL: [bool; NVERTS] = [true; NVERTS];

pub fn mark_verts(verts: &[usize]) -> [bool; NVERTS] {
    if verts.len() > NVERTS {
        panic!("verts.len() = {} > NVERTS = {}", verts.len(), NVERTS);
    }
    let mut r = [false; NVERTS];
    for &v in verts {
        r[v] = true;
    }
    r
}
