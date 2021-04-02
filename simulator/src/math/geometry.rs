// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::v2d::V2d;
use crate::math::{
    close_to_zero, in_unit_interval, max_f64s, min_f64s,
    InUnitInterval,
};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Mul};

const INTERSECTION_CLOSE_EPS: f64 = 1e-16;

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct Poly {
    pub verts: [V2d; NVERTS],
    pub edges: [LineSeg2D; NVERTS],
    pub bbox: BBox,
}

#[derive(Clone, Deserialize, Serialize)]
pub struct ExpandedPoly {
    pub verts: Vec<V2d>,
    pub edges: Vec<LineSeg2D>,
    pub bbox: BBox,
}

impl Poly {
    pub fn gen_edges(verts: &[V2d; NVERTS]) -> [LineSeg2D; NVERTS] {
        let mut edges =
            [LineSeg2D::new(&V2d::zeros(), &V2d::zeros()); NVERTS];
        (0..NVERTS).for_each(|vi| {
            edges[vi]
                .refresh(&verts[vi], &verts[circ_ix_plus(vi, NVERTS)])
        });
        edges
    }

    pub fn gen_verts(edges: &[LineSeg2D; NVERTS]) -> [V2d; NVERTS] {
        let mut verts = [V2d::default(); NVERTS];
        (0..NVERTS).for_each(|vi| {
            verts[vi] = edges[vi].p0;
        });
        verts
    }

    pub fn from_verts(verts: &[V2d; NVERTS]) -> Poly {
        let bbox = BBox::from_points(verts);
        let edges = Poly::gen_edges(verts);

        Poly {
            verts: *verts,
            edges,
            bbox,
        }
    }

    pub fn from_edges(edges: &[LineSeg2D; NVERTS]) -> Poly {
        let verts = Poly::gen_verts(edges);
        let bbox = BBox::from_points(&verts);
        Poly {
            verts,
            edges: *edges,
            bbox,
        }
    }

    pub fn expand(&self, factor: f64) -> ExpandedPoly {
        unimplemented!()
        // let mut verts = [V2d::default(); NVERTS];
        // let uevs = self
        //     .edges
        //     .iter()
        //     .map(|e| e.vector.unitize())
        //     .collect::<Vec<V2d>>();
        // let mut uivs = [V2d::default(); NVERTS];
        // (0..NVERTS).for_each(|j| {
        //     let i = circ_ix_minus(j, NVERTS);
        //     let tangent = (uevs[j] + uevs[i]).unitize();
        //     uivs[j] = tangent.normal();
        // });
        // (0..NVERTS).for_each(|i| {
        //     verts[i] = self.verts[i] + uivs[i].scale(factor);
        // });
        // Poly::from_verts(&verts)
    }
}

/// Calculate the area of a polygon with vertices positioned at `xys`.
/// [ref](http://geomalgorithms.com/a01-_area.html)
pub fn calc_poly_area(xys: &[V2d]) -> f64 {
    let nvs = xys.len();

    let mut area = 0.0_f64;
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
    pub xmin: f64,
    pub ymin: f64,
    pub xmax: f64,
    pub ymax: f64,
}

impl BBox {
    pub fn from_point_pair(a: &V2d, b: &V2d) -> BBox {
        BBox {
            xmin: a.x.min(b.x),
            ymin: a.y.min(b.y),
            xmax: a.x.max(b.x),
            ymax: a.y.max(b.y),
        }
    }

    pub fn from_points(ps: &[V2d]) -> BBox {
        let xs: Vec<f64> = ps.iter().map(|v| v.x).collect();
        let ys: Vec<f64> = ps.iter().map(|v| v.y).collect();
        BBox {
            xmin: min_f64s(&xs),
            ymin: min_f64s(&ys),
            xmax: max_f64s(&xs),
            ymax: max_f64s(&ys),
        }
    }

    #[inline]
    pub fn expand_by(&self, l: f64) -> BBox {
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
    pub fn contains(&self, point: &V2d, eps: f64) -> bool {
        point.x > self.xmin - eps
            && point.x < self.xmax + eps
            && point.y > self.ymin - eps
            && point.y < self.ymax + eps
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

pub fn is_point_in_poly(
    p: &V2d,
    poly_bbox: Option<&BBox>,
    poly: &[V2d],
) -> bool {
    if let Some(bb) = poly_bbox {
        if !bb.contains(p, INTERSECTION_CLOSE_EPS) {
            return false;
        }
    }
    let mut wn: isize = 0;
    let nverts = poly.len();
    for vi in 0..nverts {
        let p0 = &poly[vi];
        let p1 = &poly[circ_ix_plus(vi, nverts)];

        if p0.y < p.y || close_to_zero(p0.y - p.y, 1e-16) {
            if p1.y > p.y {
                if let IsLeftResult::Left =
                    is_left_pointwise(p0, p1, p)
                {
                    //println!("vi: {}, + 1", vi);
                    wn += 1;
                }
            }
        } else if p1.y < p.y || close_to_zero(p1.y - p.y, 1e-16) {
            if let IsLeftResult::Right = is_left_pointwise(p0, p1, p)
            {
                //println!("vi: {}, - 1", vi);
                wn -= 1;
            }
        }
    }
    //println!("wn: {}", wn);
    wn != 0
}

/// A line segment from p0 to p1 is the set of points `q = tp + p0`,
/// where `p = (p1 - p0)`, and `0 <= t <= 1`.
#[derive(Clone, Copy, Serialize, Deserialize, Debug, PartialEq)]
pub struct LineSeg2D {
    /// First point defining line segment.
    pub p0: V2d,
    /// Second point defining line segment.
    pub p1: V2d,
    /// `(p1 - p0)`, the vector generator of the line segment.
    pub vector: V2d,
    /// Bounding box of the segment.
    bbox: BBox,
    pub len: f64,
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

impl Add for LineSeg2D {
    type Output = LineSeg2D;

    fn add(self, rhs: Self) -> Self::Output {
        let p0 = self.p0 + rhs.p0;
        let p1 = self.p1 + rhs.p1;
        LineSeg2D::new(&p0, &p1)
    }
}

impl Mul<f64> for LineSeg2D {
    type Output = LineSeg2D;

    fn mul(self, rhs: f64) -> Self::Output {
        let p0 = self.p0.scale(rhs);
        let p1 = self.p1.scale(rhs);
        LineSeg2D::new(&p0, &p1)
    }
}

impl Mul<LineSeg2D> for f64 {
    type Output = LineSeg2D;

    fn mul(self, rhs: LineSeg2D) -> Self::Output {
        let p0 = rhs.p0.scale(self);
        let p1 = rhs.p1.scale(self);
        LineSeg2D::new(&p0, &p1)
    }
}

#[derive(Copy, Clone)]
pub enum IntersectCalcResult {
    Strict(f64, f64),
    Weak(f64, f64),
    TwoNonIdenticalPoints,
    NoBBoxOverlap,
    SelfIsPointNotOnOther,
    SelfIsCollinearPointNotOnOther,
    OtherIsCollinearPointNotOnSelf,
    OtherIsPointNotOnSelf,
    IntersectionPointNotOnSelf(f64),
    IntersectionPointNotOnOther(f64),
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

pub enum IsLeftResult {
    Left,
    Right,
    Collinear,
}

/// Consider a point `p` to be left of a line segment going from `p0` to `p1`,
/// if the cross product of the vector `p1 - p0` with `p - p0` is positive.
/// If it is negative, consider it to be to the left. If it is `0`, then it is
/// collinear with the line segment.
pub fn is_left(lseg: &LineSeg2D, p: &V2d, eps: f64) -> IsLeftResult {
    let r = p - &lseg.p0;
    let cross = lseg.vector.cross(&r);

    if close_to_zero(cross, eps) {
        IsLeftResult::Collinear
    } else if cross > 0.0 {
        IsLeftResult::Left
    } else {
        IsLeftResult::Right
    }
}

/// Version of `is_left` which takes a point-wise definition of the focus
/// line segment.
pub fn is_left_pointwise(
    p0: &V2d,
    p1: &V2d,
    p: &V2d,
) -> IsLeftResult {
    let r = p - p0;
    let vector = p1 - p0;
    let cross = vector.cross(&r);

    if close_to_zero(cross, INTERSECTION_CLOSE_EPS) {
        IsLeftResult::Collinear
    } else if cross > 0.0 {
        IsLeftResult::Left
    } else {
        IsLeftResult::Right
    }
}

pub enum CheckIntersectResult {
    Strong,
    Self0OnOther0,
    Self0OnOther1,
    Self1OnOther0,
    Self1OnOther1,
    No,
    Unknown,
}

impl LineSeg2D {
    /// Create new line segment from two given points.
    pub fn new(p0: &V2d, p1: &V2d) -> LineSeg2D {
        let p = p1 - p0;
        LineSeg2D {
            p0: *p0,
            p1: *p1,
            vector: p,
            bbox: BBox::from_point_pair(p0, p1),
            len: p.mag(),
        }
    }

    /// Refresh data in LineSeg2D points defining new line segment.
    pub fn refresh(&mut self, p0: &V2d, p1: &V2d) {
        self.p0 = *p0;
        self.p1 = *p1;
        self.vector = p1 - p0;
        self.bbox = BBox::from_point_pair(p0, p1);
    }

    /// Check if line segments intersect, without calculating point of intersection.
    pub fn check_intersection(&self, other: &LineSeg2D) -> bool {
        lsegs_intersect(&self.p0, &self.p1, other)
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
        // where `t: f64`, and `u: f64`.

        // Let `dx_s = p.x = self.p1.x - self.p0.x`; similarly, let
        // `dy_s = p.y`, `dx_o` and `dy_o` be defined analogously.
        // Also introduce `dx0_so = self.p0.x - other.p0.x` and let
        // `dy0_so` be defined analogously.

        let dx_s = self.vector.x;
        let dx_o = other.vector.x;
        let dy_s = self.vector.y;
        let dy_o = other.vector.y;

        // Let us quickly rule out some simple cases.
        let self_is_vertical =
            close_to_zero(dx_s, INTERSECTION_CLOSE_EPS);
        let self_is_horizontal =
            close_to_zero(dy_s, INTERSECTION_CLOSE_EPS);
        let other_is_vertical =
            close_to_zero(dx_o, INTERSECTION_CLOSE_EPS);
        let other_is_horizontal =
            close_to_zero(dy_o, INTERSECTION_CLOSE_EPS);
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
                if close_to_zero(ux - uy, INTERSECTION_CLOSE_EPS) {
                    match in_unit_interval(ux, INTERSECTION_CLOSE_EPS) {
                        InUnitInterval::In => {
                            IntersectCalcResult::Strict(0.0, ux)
                        }
                        InUnitInterval::One | InUnitInterval::Zero => {
                            IntersectCalcResult::Weak(0.0, ux)
                        }
                        InUnitInterval::Out => {
                            IntersectCalcResult::SelfIsCollinearPointNotOnOther
                        }
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
                if close_to_zero(tx - ty, INTERSECTION_CLOSE_EPS) {
                    match in_unit_interval(tx, INTERSECTION_CLOSE_EPS) {
                        InUnitInterval::In => {
                            IntersectCalcResult::Strict(tx, 0.0)
                        }
                        InUnitInterval::One | InUnitInterval::Zero => {
                            IntersectCalcResult::Weak(tx, 0.0)
                        }
                        InUnitInterval::Out => {
                            IntersectCalcResult::OtherIsCollinearPointNotOnSelf
                        }
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
                if close_to_zero(denominator, INTERSECTION_CLOSE_EPS)
                {
                    return IntersectCalcResult::ParallelLineSegs;
                }
                let dx0_so = self.p0.x - other.p0.x;
                let dy0_so = self.p0.y - other.p0.y;
                let t = (dy0_so * dx_o - dx0_so * dy_o) / denominator;
                let t_in_unit_interval =
                    in_unit_interval(t, INTERSECTION_CLOSE_EPS);
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
                let u_in_unit_interval =
                    in_unit_interval(u, INTERSECTION_CLOSE_EPS);
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

    #[inline]
    pub fn mag(&self) -> f64 {
        self.vector.mag()
    }
}

/// Uses the cross-product to check if a pair of points defining a line segment intersects
/// a line segment.
pub fn lsegs_intersect(
    p0: &V2d,
    p1: &V2d,
    other: &LineSeg2D,
) -> bool {
    let is_left_results = (
        is_left_pointwise(p0, p1, &other.p0),
        is_left_pointwise(p0, p1, &other.p1),
        is_left_pointwise(&other.p0, &other.p1, p0),
        is_left_pointwise(&other.p0, &other.p1, p1),
    );
    match is_left_results {
        (IsLeftResult::Left, IsLeftResult::Left, _, _)
        | (IsLeftResult::Right, IsLeftResult::Right, _, _)
        | (_, _, IsLeftResult::Left, IsLeftResult::Left)
        | (_, _, IsLeftResult::Right, IsLeftResult::Right) => false,
        (
            IsLeftResult::Collinear,
            IsLeftResult::Collinear,
            IsLeftResult::Collinear,
            IsLeftResult::Collinear,
        ) => BBox::from_point_pair(p0, p1)
            .intersects(&BBox::from_point_pair(&other.p0, &other.p1)),
        (_, _, _, _) => true,
    }
}
