// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::p2d::P2D;
use crate::math::{max_f32s, min_f32s};
use crate::utils::{circ_ix_minus, circ_ix_plus};

/// Calculate the area of a polygon with vertices positioned at `xys`. [ref](http://geomalgorithms.com/a01-_area.html)
pub fn calc_poly_area(xys: &[P2D]) -> f32 {
    let nvs = xys.len();

    let mut area = 0.0_f32;
    for i in 0..nvs {
        let j = circ_ix_plus(i, nvs);
        let k = circ_ix_minus(i, nvs);
        area += xys[i].x * (xys[j].y - xys[k].y);
    }

    area * 0.5
}

pub struct Bbox {
    pub(crate) xmin: f32,
    pub(crate) ymin: f32,
    pub(crate) xmax: f32,
    pub(crate) ymax: f32,
}

impl Bbox {
    pub fn calc(xys: &[P2D]) -> Bbox {
        let xs: Vec<f32> = xys.iter().map(|v| v.x).collect();
        let ys: Vec<f32> = xys.iter().map(|v| v.y).collect();
        Bbox {
            xmin: min_f32s(&xs),
            ymin: min_f32s(&ys),
            xmax: max_f32s(&xs),
            ymax: max_f32s(&xs),
        }
    }

    pub fn expand_by(&self, l: f32) -> Bbox {
        Bbox {
            xmin: self.xmin - l,
            ymin: self.ymin - l,
            xmax: self.xmax + l,
            ymax: self.ymax + l,
        }
    }

    pub fn intersects(&self, other: &Bbox) -> bool {
        !((self.xmin > other.xmax || self.xmax < other.xmin)
            || (self.ymin > other.ymax || self.ymax < other.ymin))
    }

    pub fn contains(&self, point: &P2D) -> bool {
        !(point.x < self.xmin || point.x > self.xmax || point.y < self.ymin || point.y > self.ymax)
    }
}

pub fn calc_line_and_seg_intersect(p0: &P2D, p1: &P2D, l0: &P2D, l: &P2D) -> f32 {
    let s = p1 - p0;
    //let t = ((l.y / l.x) * (p0.x - l0.x) + l0.y - p0.y) / (s.y - (s.x * l.y / l.x));
    (l.x * (l0.y - p0.y) - l.y * (l0.x - p0.x)) / (l.x * s.y - l.y * s.x)
}

pub fn move_inner_point_to_bdry(
    point: &P2D,
    move_vec: &P2D,
    poly_bbox: &Bbox,
    poly: &[P2D],
) -> P2D {
    if poly_bbox.contains(point) {
        let mut num_intersects: usize = 0;
        let mut intersects = [P2D::default(); 2];
        for ix in 0..poly.len() {
            let p0 = &poly[ix];
            let p1 = &poly[circ_ix_plus(ix, poly.len())];
            let t = calc_line_and_seg_intersect(p0, p1, point, move_vec);
            if (0.0 < t && t < 1.0)
                || (t < std::f32::EPSILON)
                || ((t - 1.0).abs() < std::f32::EPSILON)
            {
                intersects[num_intersects] = t * (p1 - p0) + *p0;
                num_intersects += 1;
                if num_intersects == 2 {
                    break;
                }
            }
        }

        if num_intersects == 0 || num_intersects == 2 {
            *point
        } else {
            intersects[0]
        }
    } else {
        *point
    }
}

pub fn calc_dist_point_to_seg(point: &P2D, s0: &P2D, s1: &P2D) -> f32 {
    let seg = s1 - s0;
    let rel_pt = point - s0;
    let t = (seg.x * rel_pt.x + seg.y * rel_pt.y) / (seg.x * seg.x + seg.y * seg.y);

    if (0.0 < t && t < 1.0) || (t - 1.0).abs() < f32::EPSILON || t.abs() < f32::EPSILON {
        let w = t * seg;
        let x = rel_pt - w;
        x.mag()
    } else {
        f32::INFINITY
    }
}
