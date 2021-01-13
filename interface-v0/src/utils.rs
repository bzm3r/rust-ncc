use druid::{Point, Vec2};
use rust_ncc::math::v2d::V2D;
use rust_ncc::NVERTS;

pub fn point_from(v: &V2D) -> Point {
    Point::new(v.x as f64, v.y as f64)
}

pub fn vec2_from(v: &V2D) -> Vec2 {
    Vec2::new(v.x as f64, v.y as f64)
}

pub fn point_arr_from(vs: &[V2D; NVERTS]) -> [Point; NVERTS] {
    let mut r = [Point::default(); NVERTS];
    r.iter_mut()
        .zip(vs.iter())
        .for_each(|(p, v)| *p = point_from(v));
    r
}

pub fn vec2_arr_from(vs: &[V2D; NVERTS]) -> [Vec2; NVERTS] {
    let mut r = [Vec2::default(); NVERTS];
    r.iter_mut()
        .zip(vs.iter())
        .for_each(|(p, v)| *p = vec2_from(v));
    r
}
