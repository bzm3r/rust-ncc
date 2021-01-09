// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::v2d::V2D;
use crate::math::{
    close_to_zero, in_unit_interval, max_f32s, min_f32s,
};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct Poly {
    pub verts: [V2D; NVERTS],
    pub edges: [LineSeg2D; NVERTS],
    pub bbox: BBox,
}

impl Poly {
    pub fn from_points(ps: &[V2D; NVERTS]) -> Poly {
        let bbox = BBox::from_points(ps);
        let mut edges =
            [LineSeg2D::new(&V2D::zeros(), &V2D::zeros()); NVERTS];
        (0..NVERTS).for_each(|vi| {
            edges[vi].refresh(&ps[vi], &ps[circ_ix_plus(vi, NVERTS)])
        });
        Poly {
            verts: *ps,
            edges,
            bbox,
        }
    }
}

/// Calculate the area of a polygon with vertices positioned at `xys`.
/// [ref](http://geomalgorithms.com/a01-_area.html)
pub fn calc_poly_area(xys: &[V2D]) -> f32 {
    let nvs = xys.len();

    let mut area = 0.0_f32;
    for i in 0..nvs {
        let j = circ_ix_plus(i, nvs);
        let k = circ_ix_minus(i, nvs);
        area += xys[i].x * (xys[j].y - xys[k].y);
    }

    area * 0.5
}

#[derive(Copy, Clone, Deserialize, Serialize)]
pub struct BBox {
    pub xmin: f32,
    pub ymin: f32,
    pub xmax: f32,
    pub ymax: f32,
}

impl BBox {
    pub fn from_point_pair(a: &V2D, b: &V2D) -> BBox {
        BBox {
            xmin: a.x.min(b.x),
            ymin: a.y.min(b.y),
            xmax: a.x.max(b.x),
            ymax: a.y.max(b.y),
        }
    }

    pub fn from_points(ps: &[V2D]) -> BBox {
        let xs: Vec<f32> = ps.iter().map(|v| v.x).collect();
        let ys: Vec<f32> = ps.iter().map(|v| v.y).collect();
        BBox {
            xmin: min_f32s(&xs),
            ymin: min_f32s(&ys),
            xmax: max_f32s(&xs),
            ymax: max_f32s(&ys),
        }
    }

    pub fn expand_by(&self, l: f32) -> BBox {
        BBox {
            xmin: self.xmin - l,
            ymin: self.ymin - l,
            xmax: self.xmax + l,
            ymax: self.ymax + l,
        }
    }

    pub fn intersects(&self, other: &BBox) -> bool {
        !((self.xmin > other.xmax || self.xmax < other.xmin)
            || (self.ymin > other.ymax || self.ymax < other.ymin))
    }

    pub fn contains(&self, point: &V2D) -> bool {
        !(point.x < self.xmin
            || point.x > self.xmax
            || point.y < self.ymin
            || point.y > self.ymax)
    }
}

pub enum PointSegRelation {
    Left,
    Right,
    On,
}

pub fn is_left(p: &V2D, p0: &V2D, p1: &V2D) -> PointSegRelation {
    let r =
        (p1.x - p0.x) * (p.y - p0.y) - (p.x - p0.x) * (p1.y - p0.y);
    match 0.0.partial_cmp(&r) {
        Some(Ordering::Less) => PointSegRelation::Left,
        Some(Ordering::Greater) => PointSegRelation::Right,
        Some(Ordering::Equal) => PointSegRelation::On,
        None => panic!(
            "cannot compare point {} with line seg defined ({}, {})",
            p, p0, p1
        ),
    }
}

pub fn is_point_in_poly(
    p: &V2D,
    poly_bbox: Option<&BBox>,
    poly: &[V2D],
) -> bool {
    if let Some(bb) = poly_bbox {
        if !bb.contains(p) {
            return false;
        }
    }
    let mut wn: isize = 0;
    let nverts = poly.len();
    for vi in 0..nverts {
        let p0 = &poly[vi];
        let p1 = &poly[circ_ix_plus(vi, nverts)];

        if (p0.y - p.y).abs() < f32::EPSILON || p0.y < p.y {
            if p1.y > p.y {
                if let PointSegRelation::Left = is_left(p, p0, p1) {
                    wn += 1;
                }
            }
        } else if (p1.y - p.y).abs() < f32::EPSILON || p1.y < p.y {
            if let PointSegRelation::Right = is_left(p, p0, p1) {
                wn -= 1;
            }
        }
    }
    wn != 0
}

/// A line segment from p0 to p1 is the set of points `q = tp + p0`,
/// where `p = (p1 - p0)`, and `0 <= t <= 1`.
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct LineSeg2D {
    /// First point defining line segment.
    p0: V2D,
    /// Second point defining line segment.
    p1: V2D,
    /// `(p1 - p0)`, the vector generator of the line segment.
    p: V2D,
    /// Bounding box of the segment.
    bbox: BBox,
}

impl LineSeg2D {
    /// Create new line segment from two given points.
    pub fn new(p0: &V2D, p1: &V2D) -> LineSeg2D {
        let p = p1 - p0;
        LineSeg2D {
            p0: *p0,
            p1: *p1,
            p,
            bbox: BBox::from_point_pair(p0, p1),
        }
    }

    /// Refresh data in LineSeg2D points defining new line segment.
    pub fn refresh(&mut self, p0: &V2D, p1: &V2D) {
        self.p0 = *p0;
        self.p1 = *p1;
        self.p = p1 - p0;
        self.bbox = BBox::from_point_pair(p0, p1);
    }

    // /// Create new line segment from four coordinate points
    // /// `(a, b, x, y)` of type `(f32; 4)`, assuming that `p0`
    // /// is `V2d::new(a, b)` and `p1` is `V2d::new(x, y)`.
    // pub fn from_coordinates(
    //     a: f32,
    //     b: f32,
    //     x: f32,
    //     y: f32,
    // ) -> LineSeg2D {
    //     LineSeg2D::new(&V2D::new(a, b), &V2D::new(x, y))
    // }

    /// Let this segment `self` be parametrized so that a point
    /// `p` lies on `self` if `p = t * self.p + self.p0` for
    /// `0 <= t <= 1`. Similarly, the `other` line segment be
    /// parametrized so that a point `q` lies on `other` if  
    /// `p = u * self.p + self.p0` for `0 <= u <= 1`.
    ///
    /// This function calculates `(t, u)`, if an intersection exists.
    pub fn intersects_lseg(
        &self,
        other: &LineSeg2D,
    ) -> Option<(f32, f32)> {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return None;
        }
        // Let `s` be a point on this ("self") line segment, and let
        // `o` be a point on the "other" line segment, such that:
        // `s = t * self.p + self.p0`
        // `o = u * self.p + self.p0`
        // where `t: f32`, and `u: f32`.

        // Let `dx_s = p.x = self.p1.x - self.p0.x`; similarly, let
        // `dy_s = p.y`, `dx_o` and `dy_o` be defined analogously.
        // Also introduce `dx0_so = self.p0.x - other.p0.x` and let
        // `dy0_so` be defined analgoously.

        let dx_s = self.p.x;
        let dx_o = other.p.x;
        let dy_s = self.p.y;
        let dy_o = other.p.y;

        // Let us quickly rule out some simple cases.
        let self_is_vertical = close_to_zero(dx_s);
        let self_is_horizontal = close_to_zero(dy_s);
        let other_is_vertical = close_to_zero(dx_o);
        let other_is_horizontal = close_to_zero(dy_o);
        match (
            self_is_vertical,
            self_is_horizontal,
            other_is_vertical,
            other_is_horizontal,
        ) {
            (true, true, true, true) => {
                // Both line segments are generated by a zero vector. In other words,
                // they are indistinguishable from points. So, just check if points
                // are the same.
                if (self.p0 - self.p1).close_to_zero() {
                    Some((0.0, 0.0))
                } else {
                    None
                }
            }
            (true, true, _, _) => {
                // `self` is a point, so just check if it lies
                // `other`.
                let ux = (self.p0.x - other.p0.x) / dx_o;
                let uy = (self.p0.y - other.p0.y) / dy_o;
                if in_unit_interval(ux)
                    && in_unit_interval(uy)
                    && close_to_zero(ux - uy)
                {
                    Some((0.0, ux))
                } else {
                    None
                }
            }
            (_, _, true, true) => {
                // `other` is a point, so just check if it lies
                // `self`.
                let tx = (other.p0.x - self.p0.x) / dx_o;
                let ty = (other.p0.y - self.p0.y) / dy_o;
                if in_unit_interval(tx)
                    && in_unit_interval(ty)
                    && close_to_zero(tx - ty)
                {
                    Some((0.0, tx))
                } else {
                    None
                }
            }
            (_, _, _, _) => {
                // Note that `s.x = t * dx_s + self.p0.x`
                //           `s.y = t * dy_s + self.p0.y`
                // Similarly `o.x = u * dx_o + other.p0.x`
                //           `o.y = u * dy_o + other.p0.y`
                // At intersection, we have `s = o`, that is, `s.x = o.x` and
                // `s.y = o.y`. So:
                //      `t * dx_s + self.p0.x = u * dx_o + other.p0.x`
                //      { `let dx0_so = self.p0.x - other.p0.x` }
                // <=>  `t * (dx_s / dx_o) + (dx0_so / dx_o) =  u`
                // <=>  `t * (dx_s / dx_o) + (dx0_so / dx_o) =  u`
                // And similarly we have:
                //      `t * (dy_s / dy_o) + (dy0_so / dy_o) = u`
                // Therefore:
                //      `t * (dx_s / dx_o) + (dx0_so / dx_o) = t * (dy_s / dy_o) + (dy0_so / dy_o)`
                //      { re-arrange terms with `t` onto one side, and terms without on the other }
                // <=>  `t * (dx_s / dx_o) - t * (dy_s / dy_o) = (dy0_so / dy_o) - (dx0_so / dx_o)`
                //      { factor out `t`, add the fractions on both sides }
                // <=>  `t * (dx_s * dy_o - dy_s * dx_o)/ (dx_o * dy_o) = (dy0_so * dx_o - dx0_so * dy_o)/ (dx_o * dy_o)`
                //      { cancel out `(dx_o * dy_o)` on both sides, and isolate `t` }
                // <=>  `t = (dy0_so * dx_o - dx0_so * dy_o)/(dx_s * dy_o - dy_s * dx_o)`
                //
                // Notice that the denominator is the cross-product of the
                // two vectors. This product is zero precisely when generating
                // vectors of both line segments are parallel. Therefore,
                // let us check to see if that is the case.
                let denominator = dx_s * dy_o - dy_s * dx_o;
                if close_to_zero(denominator) {
                    return None;
                }
                let dx0_so = self.p0.x - other.p0.x;
                let dy0_so = self.p0.y - other.p0.y;
                let t = (dy0_so * dx_o - dx0_so * dy_o) / denominator;
                if !in_unit_interval(t) {
                    return None;
                }
                // Recalling that `u = t * (dx_s / dx_o) + (dx0_so / dx_o)`
                // and `t * (dy_s / dy_o) + (dy0_so / dy_o)`. Which formula
                // we use to determine `u` depends on whether or not `other`
                // is vertical.
                let u = match (other_is_vertical, other_is_horizontal)
                {
                    (false, _) => t * (dx_s / dx_o) + (dx0_so / dx_o),
                    (_, false) => t * (dy_s / dy_o) + (dy0_so / dy_o),
                    (true, true) => {
                        panic!(
                            "Did not expect to reach case where\
                         `other` is a point. It should have been\
                          handled by earlier match cases"
                        );
                    }
                };
                if !in_unit_interval(u) {
                    return None;
                }
                Some((t, u))
            }
        }
    }

    // pub fn calc_intersect_point(
    //     &self,
    //     p0: &V2D,
    //     p1: &V2D,
    // ) -> Option<V2D> {
    //     let other = LineSeg2D::new(p0, p1);
    //     self.intersects_lseg(&other)
    //         .map(|(t, _)| t * self.p + self.p0)
    // }

    pub fn intersects_bbox(&self, bbox: &BBox) -> bool {
        self.bbox.intersects(bbox)
    }

    pub fn intersects_poly(&self, poly: &Poly) -> bool {
        if !self.intersects_bbox(&poly.bbox) {
            return false;
        }

        for edge in poly.edges.iter() {
            match self.intersects_lseg(&edge) {
                Some((t, u)) => {
                    match (t > 0.0, t < 1.0, u > 0.0, u < 1.0) {
                        (true, true, true, true) => {
                            return true;
                        }
                        _ => {
                            continue;
                        }
                    }
                }
                None => {
                    continue;
                }
            }
        }
        false
    }

    #[inline]
    pub fn mag(&self) -> f32 {
        self.p.mag()
    }
}

// pub fn refine_raw_poly(raw_poly: [[f32; 2]; 16]) -> [V2d; 16] {
//     let mut r = [V2d::default(); 16];
//     for (q, p) in r.iter_mut().zip(raw_poly.iter()) {
//         q.x = p[0];
//         q.y = p[1];
//     }
//     r
// }
//
// pub fn are_polys_equal(poly0: &[V2d; 16], poly1: &[V2d; 16]) -> bool {
//     for (p, q) in poly0.iter().zip(poly1.iter()) {
//         if p.x != q.x || p.y != q.y {
//             return false;
//         }
//     }
//     true
// }
//
// pub fn debug_point_in_poly() {
//     let c0_347 = refine_raw_poly([
//         [40.927586, 17.505148],
//         [38.242447, 24.997238],
//         [32.53913, 30.540205],
//         [26.321924, 35.535126],
//         [19.893433, 40.224964],
//         [15.902101, 36.503548],
//         [9.020461, 31.977158],
//         [3.0782228, 26.634775],
//         [-0.060920116, 19.310455],
//         [0.8173119, 11.377405],
//         [5.091162, 4.621307],
//         [11.693666, 0.09866802],
//         [19.504337, -1.6616803],
//         [27.426825, -0.5139225],
//         [34.400642, 3.4034972],
//         [39.32033, 9.69614],
//     ]);
//     let c1_348: [V2d; 16] = refine_raw_poly([
//         [40.00125, 61.21729],
//         [39.038857, 69.13645],
//         [34.721806, 75.863235],
//         [28.110025, 80.3738],
//         [20.295544, 82.1296],
//         [12.370754, 80.97098],
//         [5.389173, 77.058975],
//         [0.44813222, 70.77858],
//         [-1.1411942, 62.963577],
//         [1.6184937, 55.497562],
//         [7.390409, 50.025238],
//         [13.649933, 45.080173],
//         [17.333612, 37.841843],
//         [24.296999, 43.844837],
//         [31.142982, 48.41107],
//         [36.987385, 53.847706],
//     ]);
//     let c1_348_in_c0_check: [V2d; 16] = refine_raw_poly([
//         [40.00125, 61.21729],
//         [39.038857, 69.13645],
//         [34.721806, 75.863235],
//         [28.110025, 80.3738],
//         [20.295544, 82.1296],
//         [12.370754, 80.97098],
//         [5.389173, 77.058975],
//         [0.44813222, 70.77858],
//         [-1.1411942, 62.963577],
//         [1.6184937, 55.497562],
//         [7.390409, 50.025238],
//         [13.649933, 45.080173],
//         [17.333612, 37.841843],
//         [24.296999, 43.844837],
//         [31.142982, 48.41107],
//         [36.987385, 53.847706],
//     ]);
//     let c1_348_is_ok = are_polys_equal(&c1_348, &c1_348_in_c0_check);
//     let any_c1_348_in_c0_347 = c1_348
//         .iter()
//         .any(|p| is_point_in_poly_no_bb_check(p, &c0_347));
//     let c0v4_in_c1 =
//         is_point_in_poly_no_bb_check(&c0_347[4], &c1_348);
//     let c0v4_in_c1ve =
//         is_point_in_poly_no_bb_check(&c0_347[4], &c1_348_in_c0_check);
//     println!("is c1_348 okay: {}", c1_348_is_ok);
//     println!(
//         "any vertex of c1_348 in c0_347: {}",
//         any_c1_348_in_c0_347
//     );
//     println!("is c0_347_v4 in c1_348: {}", c0v4_in_c1);
// }
