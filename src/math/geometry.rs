// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::v2d::V2d;
use crate::math::{max_f32s, min_f32s};
use crate::utils::{
    circ_ix_minus, circ_ix_plus,
};
use std::cmp::Ordering;

/// Calculate the area of a polygon with vertices positioned at `xys`. [ref](http://geomalgorithms.com/a01-_area.html)
pub fn calc_poly_area(
    xys: &[V2d],
) -> f32 {
    let nvs = xys.len();

    let mut area = 0.0_f32;
    for i in 0..nvs {
        let j = circ_ix_plus(i, nvs);
        let k = circ_ix_minus(i, nvs);
        area += xys[i].x
            * (xys[j].y - xys[k].y);
    }

    area * 0.5
}

#[derive(Copy, Clone)]
pub struct BBox {
    pub xmin: f32,
    pub ymin: f32,
    pub xmax: f32,
    pub ymax: f32,
}

impl BBox {
    pub fn from_points(
        xys: &[V2d],
    ) -> BBox {
        let xs: Vec<f32> = xys
            .iter()
            .map(|v| v.x)
            .collect();
        let ys: Vec<f32> = xys
            .iter()
            .map(|v| v.y)
            .collect();
        BBox {
            xmin: min_f32s(&xs),
            ymin: min_f32s(&ys),
            xmax: max_f32s(&xs),
            ymax: max_f32s(&ys),
        }
    }

    pub fn expand_by(
        &self,
        l: f32,
    ) -> BBox {
        BBox {
            xmin: self.xmin - l,
            ymin: self.ymin - l,
            xmax: self.xmax + l,
            ymax: self.ymax + l,
        }
    }

    pub fn intersects(
        &self,
        other: &BBox,
    ) -> bool {
        !((self.xmin > other.xmax
            || self.xmax < other.xmin)
            || (self.ymin > other.ymax
                || self.ymax
                    < other.ymin))
    }

    pub fn contains(
        &self,
        point: &V2d,
    ) -> bool {
        !(point.x < self.xmin
            || point.x > self.xmax
            || point.y < self.ymin
            || point.y > self.ymax)
    }
}

#[allow(unused)]
fn in_unit_interval(x: f32) -> bool {
    (x.abs() < f32::EPSILON
        || (x - 1.0).abs()
            < f32::EPSILON)
        || (0.0 < x && x < 1.0)
}

pub enum PointSegRelation {
    Left,
    Right,
    On,
}

pub fn is_left(
    p: &V2d,
    p0: &V2d,
    p1: &V2d,
) -> PointSegRelation {
    let r = (p1.x - p0.x)
        * (p.y - p0.y)
        - (p.x - p0.x) * (p1.y - p0.y);
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
    p: &V2d,
    poly_bbox: &BBox,
    poly: &[V2d],
) -> bool {
    if poly_bbox.contains(p) {
        is_point_in_poly_no_bb_check(
            p, poly,
        )
    } else {
        false
    }
}

pub fn is_point_in_poly_no_bb_check(
    p: &V2d,
    poly: &[V2d],
) -> bool {
    let mut wn: isize = 0;
    let nverts = poly.len();
    for vi in 0..nverts {
        let p0 = &poly[vi];
        let p1 = &poly
            [circ_ix_plus(vi, nverts)];

        if (p0.y - p.y).abs()
            < f32::EPSILON
            || p0.y < p.y
        {
            if p1.y > p.y {
                if let PointSegRelation::Left = is_left(p, p0, p1) {
                    wn += 1;
                }
            }
        } else if (p1.y - p.y).abs()
            < f32::EPSILON
            || p1.y < p.y
        {
            if let PointSegRelation::Right = is_left(p, p0, p1) {
                wn -= 1;
            }
        }
    }
    wn != 0
}

/// Returns (t, d), where `k = (s1 - s0)*t + s1` is the point on `s0` to `s1` closest to `point`.
pub fn calc_dist_point_to_seg(
    point: &V2d,
    s0: &V2d,
    s1: &V2d,
) -> (f32, f32) {
    let seg = s1 - s0;
    let rel_pt = point - s0;
    let t = (seg.x * rel_pt.x
        + seg.y * rel_pt.y)
        / (seg.x * seg.x
            + seg.y * seg.y);

    if t < 0.0 || t > 1.0 {
        (f32::INFINITY, f32::INFINITY)
    } else {
        let w = t * seg;
        let x = rel_pt - w;
        (t, x.mag())
    }
}

pub struct LineSeg {
    p0: V2d,
    p1: V2d,
    p: V2d,
}

impl LineSeg {
    pub fn new(
        p0: &V2d,
        p1: &V2d,
    ) -> LineSeg {
        let p = p1 - p0;
        LineSeg {
            p0: *p0,
            p1: *p1,
            p,
        }
    }

    pub fn check_intersection(
        &self,
        p0: &V2d,
        p1: &V2d,
    ) -> Option<f32> {
        // TODO: make sure this is correct!
        let p = p0 - p1;
        let alpha = self.p0.x - p0.x;
        let beta = self.p1.y - p1.y;
        let gamma = self.p.y * p.x
            - self.p.x * p.y;
        let t = -1.0
            * (p.y * alpha
                - p.x * beta)
            / gamma;
        let u = self.p.y * alpha
            - self.p.x * beta / gamma;

        match in_unit_interval(t)
            && in_unit_interval(u)
        {
            true => Some(t),
            false => None,
        }
    }

    #[allow(unused)]
    pub fn calc_intersection(
        &self,
        p0: &V2d,
        p1: &V2d,
    ) -> Option<V2d> {
        match self
            .check_intersection(p0, p1)
        {
            Some(t) => Some(
                self.p0 + t * self.p,
            ),
            None => None,
        }
    }
}

pub fn ls_self_intersects_poly(
    vi: usize,
    poly: &[V2d],
    lseg: &LineSeg,
) -> bool {
    let nverts = poly.len();
    poly.iter().enumerate().any(
        |(i, vc)| {
            let j =
                circ_ix_plus(i, nverts);
            if i != vi && j != vi {
                let vc1 = &poly[j];
                lseg.check_intersection(
                    vc, vc1,
                )
                .is_some()
            } else {
                false
            }
        },
    )
}

pub fn ls_intersects_poly(
    lseg: &LineSeg,
    poly: &[V2d],
) -> bool {
    let nverts = poly.len();
    poly.iter().enumerate().any(
        |(vi, vc)| {
            let vc1 = &poly
                [circ_ix_plus(
                    vi, nverts,
                )];
            lseg.check_intersection(
                vc, vc1,
            )
            .is_some()
        },
    )
}
