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
    InUnitInterval,
};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt;
use std::fmt::Display;

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

#[derive(
    Copy, Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
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

    #[inline]
    pub fn expand_by(&self, l: f32) -> BBox {
        BBox {
            xmin: self.xmin - l,
            ymin: self.ymin - l,
            xmax: self.xmax + l,
            ymax: self.ymax + l,
        }
    }

    #[inline]
    pub fn intersects(&self, other: &BBox) -> bool {
        self.xmax > other.xmin
            && self.xmin < other.xmax
            && self.ymin < other.ymax
            && self.ymax > other.ymin
    }

    #[inline]
    pub fn contains(&self, point: &V2D) -> bool {
        point.x > self.xmin
            && point.x < self.xmax
            && point.y > self.ymin
            && point.y < self.ymax
    }
}

impl Display for BBox {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "BBox {{ bot_left: ({}, {}), top_right: ({}, {}) }}",
            self.xmin, self.ymin, self.xmax, self.ymax
        )
    }
}

pub enum PointSegRelation {
    Left,
    Right,
    Collinear,
}

pub fn is_left(p: &V2D, p0: &V2D, p1: &V2D) -> PointSegRelation {
    let r =
        (p1.x - p0.x) * (p.y - p0.y) - (p.x - p0.x) * (p1.y - p0.y);
    match 0.0.partial_cmp(&r) {
        Some(Ordering::Less) => PointSegRelation::Left,
        Some(Ordering::Greater) => PointSegRelation::Right,
        Some(Ordering::Equal) => PointSegRelation::Collinear,
        None => panic!(
            "cannot compare point {} with line seg defined ({}, {})",
            p, p0, p1
        ),
    }
}

#[inline]
pub fn is_left_simple(p: &V2D, p0: &V2D, p1: &V2D) -> i32 {
    let r =
        (p1.x - p0.x) * (p.y - p0.y) - (p.x - p0.x) * (p1.y - p0.y);
    if r > 0.0 {
        1
    } else if r < 0.0 {
        -1
    } else {
        0
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
    pub vector: V2D,
    /// Bounding box of the segment.
    bbox: BBox,
}

impl Display for LineSeg2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "LineSeg2D {{ p0: {}, p1: {}, bbox: {} }}",
            self.p0, self.p1, self.bbox
        )
    }
}

#[derive(Copy, Clone)]
pub enum IntersectCalcResult {
    Strict(f32, f32),
    Weak(f32, f32),
    TwoNonIdenticalPoints,
    NoBBoxOverlap,
    SelfIsPointNotOnOther,
    SelfIsCollinearPointNotOnOther,
    OtherIsCollinearPointNotOnSelf,
    OtherIsPointNotOnSelf,
    IntersectionPointNotOnSelf(f32),
    IntersectionPointNotOnOther(f32),
    ParallelLineSegs,
}

impl Display for IntersectCalcResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match &self {
            IntersectCalcResult::Strict(t, u) => {
                format!("Intersection({}, {})", t, u)
            }
            IntersectCalcResult::Weak(t, u) => {
                format!("WeakIntersection({}, {})", t, u)
            }
            IntersectCalcResult::IntersectionPointNotOnSelf(t) => {
                format!("IntersectionPointNotOnSelf({})", t)
            }
            IntersectCalcResult::IntersectionPointNotOnOther(u) => {
                format!("IntersectionPointNotOnSelf({})", u)
            }
            IntersectCalcResult::NoBBoxOverlap => {
                "NoBBoxOverlap".to_string()
            }
            IntersectCalcResult::OtherIsPointNotOnSelf => {
                "OtherIsPointNotOnSelf".to_string()
            }
            IntersectCalcResult::ParallelLineSegs => {
                "ParallelLineSegs".to_string()
            }
            IntersectCalcResult::SelfIsPointNotOnOther => {
                "SelfIsPointNotOnOther".to_string()
            }
            IntersectCalcResult::TwoNonIdenticalPoints => {
                "TwoNonIdenticalPoints".to_string()
            }
            IntersectCalcResult::SelfIsCollinearPointNotOnOther => {
                "SelfIsCollinearPointNotOnOther".to_string()
            }
            IntersectCalcResult::OtherIsCollinearPointNotOnSelf => {
                "OtherIsCollinearPointNotOnSelf".to_string()
            }
        };
        write!(f, "{}", s)
    }
}

impl LineSeg2D {
    /// Create new line segment from two given points.
    pub fn new(p0: &V2D, p1: &V2D) -> LineSeg2D {
        let p = p1 - p0;
        LineSeg2D {
            p0: *p0,
            p1: *p1,
            vector: p,
            bbox: BBox::from_point_pair(p0, p1),
        }
    }

    /// Refresh data in LineSeg2D points defining new line segment.
    pub fn refresh(&mut self, p0: &V2D, p1: &V2D) {
        self.p0 = *p0;
        self.p1 = *p1;
        self.vector = p1 - p0;
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

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return false;
        }

        match (
            is_left(&self.p0, &self.p1, &other.p0),
            is_left(&self.p0, &self.p1, &other.p1),
        ) {
            (PointSegRelation::Collinear, _)
            | (_, PointSegRelation::Collinear)
            | (PointSegRelation::Left, PointSegRelation::Left)
            | (PointSegRelation::Right, PointSegRelation::Right) => {
                return false;
            }
            _ => {}
        }

        match (
            is_left(&other.p0, &other.p1, &self.p0),
            is_left(&other.p0, &other.p1, &self.p1),
        ) {
            (PointSegRelation::Collinear, _)
            | (_, PointSegRelation::Collinear)
            | (PointSegRelation::Left, PointSegRelation::Left)
            | (PointSegRelation::Right, PointSegRelation::Right) => {
                return false;
            }
            _ => {}
        }

        true
    }

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect_v2(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return false;
        }

        let is_left_op0 =
            is_left_simple(&self.p0, &self.p1, &other.p0);
        if is_left_op0 == 0 {
            return false;
        }
        let is_left_op1 =
            is_left_simple(&self.p0, &self.p1, &other.p1);
        if is_left_op1 == 0 {
            return false;
        }
        if is_left_op0 == is_left_op1 {
            return false;
        }

        let is_left_sp0 =
            is_left_simple(&other.p0, &other.p1, &self.p0);
        if is_left_sp0 == 0 {
            return false;
        }
        let is_left_sp1 =
            is_left_simple(&other.p0, &other.p1, &self.p1);
        if is_left_sp1 == 0 {
            return false;
        }
        if is_left_sp0 == is_left_sp1 {
            return false;
        }

        true
    }

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect_v3(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return false;
        }
        let p0 = self.p0;
        let p1 = self.p1;
        let q0 = other.p0;
        let q1 = other.p1;

        let is_left_q0 = (p1.x - p0.x) * (q0.y - p0.y)
            - (q0.x - p0.x) * (p1.y - p0.y);
        if is_left_q0.abs() < 1e-4 {
            return false;
        }
        let is_left_q1 = (p1.x - p0.x) * (q1.y - p0.y)
            - (q1.x - p0.x) * (p1.y - p0.y);
        if is_left_q1.abs() < 1e-4
            || (is_left_q0 - is_left_q1).abs() < 1e-4
        {
            return false;
        }
        let is_left_p0 = (q1.x - q0.x) * (p0.y - q0.y)
            - (p0.x - q0.x) * (q1.y - q0.y);
        if is_left_p0.abs() < 1e-4 {
            return false;
        }
        let is_left_p1 = (q1.x - q0.x) * (p1.y - q0.y)
            - (p1.x - q0.x) * (q1.y - q0.y);
        if is_left_p1.abs() < 1e-4
            || (is_left_p0 - is_left_p1).abs() < 1e-4
        {
            return false;
        }
        true
    }

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect_v4(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return false;
        }
        let p0 = self.p0;
        let p1 = self.p1;
        let q0 = other.p0;
        let q1 = other.p1;

        let delta_px = p1.x - p0.x;
        let delta_py = p1.y - p0.y;
        let delta0y = q0.y - p0.y;
        let delta0x = q0.x - p0.x;
        let is_left_q0 = delta_px * delta0y - delta0x * delta_py;
        if is_left_q0.abs() < 1e-4 {
            return false;
        }
        let is_left_q1 =
            delta_px * (q1.y - p0.y) - (q1.x - p0.x) * delta_py;
        if is_left_q1.abs() < 1e-4
            || (is_left_q0 - is_left_q1).abs() < 1e-4
        {
            return false;
        }
        let delta_qx = q1.x - q0.x;
        let delta_qy = q1.y - q0.y;
        let is_left_p0 = delta0x * delta_qy - delta_qx * delta0y;
        if is_left_p0.abs() < 1e-4 {
            return false;
        }
        let is_left_p1 =
            delta_qx * (p1.y - q0.y) - (p1.x - q0.x) * delta_qy;
        if is_left_p1.abs() < 1e-4
            || (is_left_p0 - is_left_p1).abs() < 1e-4
        {
            return false;
        }
        true
    }

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect_v5(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        // First check to make sure that the bounding boxes intersect.
        if self.bbox.xmin > other.bbox.xmax
            || self.bbox.xmax < other.bbox.xmin
            || self.bbox.ymin > other.bbox.ymax
            || self.bbox.ymax < other.bbox.ymin
        {
            return false;
        }

        let p0 = self.p0;
        let p1 = self.p1;
        let q0 = other.p0;

        let delta_px = p1.x - p0.x;
        let delta_py = p1.y - p0.y;
        let delta0y = q0.y - p0.y;
        let delta0x = q0.x - p0.x;

        let is_left_q0 = delta_px * delta0y - delta0x * delta_py;
        if is_left_q0.abs() < 1e-4 {
            return false;
        }

        let q1 = other.p1;
        let is_left_q1 =
            delta_px * (q1.y - p0.y) - (q1.x - p0.x) * delta_py;
        if is_left_q1.abs() < 1e-4
            || (is_left_q0 - is_left_q1).abs() < 1e-4
        {
            return false;
        }

        let delta_qx = q1.x - q0.x;
        let delta_qy = q1.y - q0.y;
        let is_left_p0 = delta0x * delta_qy - delta_qx * delta0y;
        if is_left_p0.abs() < 1e-4 {
            return false;
        }

        let is_left_p1 =
            delta_qx * (p1.y - q0.y) - (p1.x - q0.x) * delta_qy;
        if is_left_p1.abs() < 1e-4
            || (is_left_p0 - is_left_p1).abs() < 1e-4
        {
            return false;
        }

        true
    }

    /// Check to see if an intersection exists between this segment
    /// and another. Does not calculate the exact location of the
    /// intersection; for that, use `LineSeg2D::calc_lseg_intersect`.
    pub fn check_strict_lseg_intersect_v6(
        &self,
        other: &LineSeg2D,
    ) -> bool {
        let p0 = self.p0;
        let p1 = self.p1;
        let q0 = other.p0;
        let q1 = other.p1;

        let delta_px = p1.x - p0.x;
        let delta_py = p1.y - p0.y;
        let delta0y = q0.y - p0.y;
        let delta0x = q0.x - p0.x;
        let is_left_q0 = delta_px * delta0y - delta0x * delta_py;
        if is_left_q0.abs() < 1e-4 {
            return false;
        }
        let is_left_q1 =
            delta_px * (q1.y - p0.y) - (q1.x - p0.x) * delta_py;
        if is_left_q1.abs() < 1e-4
            || (is_left_q0 - is_left_q1).abs() < 1e-4
        {
            return false;
        }
        let delta_qx = q1.x - q0.x;
        let delta_qy = q1.y - q0.y;
        let is_left_p0 = delta0x * delta_qy - delta_qx * delta0y;
        if is_left_p0.abs() < 1e-4 {
            return false;
        }
        let is_left_p1 =
            delta_qx * (p1.y - q0.y) - (p1.x - q0.x) * delta_qy;
        if is_left_p1.abs() < 1e-4
            || (is_left_p0 - is_left_p1).abs() < 1e-4
        {
            return false;
        }
        true
    }

    /// Let this segment `self` be parametrized so that a point
    /// `p` lies on `self` if `p = t * self.p + self.p0` for
    /// `0 <= t <= 1`. Similarly, the `other` line segment be
    /// parametrized so that a point `q` lies on `other` if  
    /// `p = u * self.p + self.p0` for `0 <= u <= 1`.
    ///
    /// This function calculates `(t, u)`, if an intersection exists.
    pub fn calc_lseg_intersect(
        &self,
        other: &LineSeg2D,
    ) -> IntersectCalcResult {
        // First check to make sure that the bounding boxes intersect.
        if !self.bbox.intersects(&other.bbox) {
            return IntersectCalcResult::NoBBoxOverlap;
        }
        // Let `s` be a point on this ("self") line segment, and let
        // `o` be a point on the "other" line segment, such that:
        // `s = t * self.p + self.p0`
        // `o = u * self.p + self.p0`
        // where `t: f32`, and `u: f32`.

        // Let `dx_s = p.x = self.p1.x - self.p0.x`; similarly, let
        // `dy_s = p.y`, `dx_o` and `dy_o` be defined analogously.
        // Also introduce `dx0_so = self.p0.x - other.p0.x` and let
        // `dy0_so` be defined analogously.

        let dx_s = self.vector.x;
        let dx_o = other.vector.x;
        let dy_s = self.vector.y;
        let dy_o = other.vector.y;

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
                    IntersectCalcResult::Strict(0.0, 0.0)
                } else {
                    IntersectCalcResult::TwoNonIdenticalPoints
                }
            }
            (true, true, _, _) => {
                // `self` is a point, so just check if it lies on other
                // `other`.
                let ux = (self.p0.x - other.p0.x) / dx_o;
                let uy = (self.p0.y - other.p0.y) / dy_o;
                if close_to_zero(ux - uy) {
                    match in_unit_interval(ux) {
                        InUnitInterval::In => IntersectCalcResult::Strict(0.0, ux),
                        InUnitInterval::One | InUnitInterval::Zero => IntersectCalcResult::Weak(0.0, ux),
                        InUnitInterval::Out => IntersectCalcResult::SelfIsCollinearPointNotOnOther,
                    }
                } else {
                    IntersectCalcResult::SelfIsPointNotOnOther
                }
            }
            (_, _, true, true) => {
                // `other` is a point, so just check if it lies
                // `self`.
                let tx = (other.p0.x - self.p0.x) / dx_o;
                let ty = (other.p0.y - self.p0.y) / dy_o;
                if close_to_zero(tx - ty) {
                    match in_unit_interval(tx) {
                        InUnitInterval::In => IntersectCalcResult::Strict(tx, 0.0),
                        InUnitInterval::One | InUnitInterval::Zero => IntersectCalcResult::Weak(tx, 0.0),
                        InUnitInterval::Out => IntersectCalcResult::OtherIsCollinearPointNotOnSelf,
                    }
                } else {
                    IntersectCalcResult::OtherIsPointNotOnSelf
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
                    return IntersectCalcResult::ParallelLineSegs;
                }
                let dx0_so = self.p0.x - other.p0.x;
                let dy0_so = self.p0.y - other.p0.y;
                let t = (dy0_so * dx_o - dx0_so * dy_o) / denominator;
                let t_in_unit_interval = in_unit_interval(t);
                if let InUnitInterval::Out = t_in_unit_interval {
                    return IntersectCalcResult::IntersectionPointNotOnSelf(t);
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
                            "Unexpected: reached branch where other is point. Should have been handled earlier."
                        );
                    }
                };
                let u_in_unit_interval = in_unit_interval(u);
                if let InUnitInterval::Out = u_in_unit_interval {
                    return IntersectCalcResult::IntersectionPointNotOnOther(u);
                }
                match (t_in_unit_interval, u_in_unit_interval) {
                    (InUnitInterval::Zero, _)
                    | (InUnitInterval::One, _)
                    | (_, InUnitInterval::Zero)
                    | (_, InUnitInterval::One) => {
                        IntersectCalcResult::Weak(t, u)
                    }
                    (InUnitInterval::In, InUnitInterval::In) => {
                        IntersectCalcResult::Strict(t, u)
                    },
                    (InUnitInterval::Out, _) => panic!("Reached case where t not in unit interval! Should have been handled earlier."),
                    (_, InUnitInterval::Out) => panic!("Reached case where u not in unit interval! Should have been handled earlier."),
                }
            }
        }
    }

    pub fn intersects_bbox(&self, bbox: &BBox) -> bool {
        self.bbox.intersects(bbox)
    }

    pub fn check_poly_intersect(&self, poly: &Poly) -> bool {
        if !self.intersects_bbox(&poly.bbox) {
            return false;
        }

        for edge in poly.edges.iter() {
            if self.check_strict_lseg_intersect_v6(&edge) {
                return true;
            }
        }
        false
    }

    pub fn old_check_poly_intersect(&self, poly: &Poly) -> bool {
        if !self.intersects_bbox(&poly.bbox) {
            return false;
        }

        for edge in poly.edges.iter() {
            match self.calc_lseg_intersect(&edge) {
                IntersectCalcResult::Strict(_, _) => {
                    return true;
                }
                _ => {
                    continue;
                }
            }
        }
        false
    }

    #[inline]
    pub fn mag(&self) -> f32 {
        self.vector.mag()
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
