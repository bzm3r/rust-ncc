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

/// Calculate the area of a polgon with vertices positioned at `xys`. [ref](http://geomalgorithms.com/a01-_area.html)
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
}

pub fn calc_line_and_seg_intersect(p0: &P2D, p1: &P2D, l0: &P2D, l: &P2D) -> f32 {
    let s = p1 - p0;
    //let t = ((l.y / l.x) * (p0.x - l0.x) + l0.y - p0.y) / (s.y - (s.x * l.y / l.x));
    (l.x * (l0.y - p0.y) - l.y * (l0.x - p0.x)) / (l.x * s.y - l.y * s.x)
}

pub fn move_inner_point_to_bdry(point: &P2D, move_vec: &P2D, pol: &[P2D]) -> P2D {
    let mut num_intersects: usize = 0;
    let mut intersects = [P2D::default(); 2];
    for ix in 0..pol.len() {
        let p0 = &pol[ix];
        let p1 = &pol[circ_ix_plus(ix, pol.len())];
        let t = calc_line_and_seg_intersect(p0, p1, point, move_vec);
        if (0.0 < t && t < 1.0) || (t < std::f32::EPSILON) || ((t - 1.0).abs() < std::f32::EPSILON)
        {
            intersects[num_intersects] = t * (p1 - p0) + p0;
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
}
