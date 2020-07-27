use crate::utils::{circ_ix_minus, circ_ix_plus};
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;
use std::ops::{Add, Div, Mul, Sub};
use std::fmt::Display;
use std::fmt;

/// Value always between `[0.0, 2*PI]`.
#[derive(PartialOrd, PartialEq, Clone, Copy)]
pub struct Radians(f32);

const RAD_2PI: Radians = Radians(2.0 * PI);
const RAD_EPS: Radians = Radians(1e-12);

/// Addition modulo `2*PI`.
impl Add for Radians {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Radians(modulo_f32(self.0 + other.0, RAD_2PI.0))
    }
}

/// Subtraction modulo `2*PI`.
impl Sub for Radians {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Radians(modulo_f32(self.0 - other.0, RAD_2PI.0))
    }
}

impl Radians {
    /// Helper function for `Radians::between`.
    fn _between(&self, t0: Radians, t1: Radians) -> bool {
        t0 < *self && *self < t1
    }

    /// Calculates if `self` is between `t0` and `t1`, where `t1` is assumed to be counter-clockwise after `t0`.
    pub fn between(&self, t0: Radians, t1: Radians) -> bool {
        if t0 - t1 < RAD_EPS {
            false
        } else if t1 - *self < RAD_EPS || t0 - *self < RAD_EPS {
            true
        } else if t0 < t1 {
            self._between(t0, t1)
        } else {
            !self._between(t1, t0)
        }
    }
}

/// 2D vector with `f32` elements.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct P2D {
    pub x: f32,
    pub y: f32,
}

impl P2D {
    /// Calculate direction of vector from origin to point in radians.
    pub fn direction(&self) -> Radians {
        arctan(self.x, self.y)
    }

    /// Calculate magnitude of vector from origin to point.
    pub fn mag(&self) -> f32 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    pub fn dot(&self, other: &P2D) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn abs(&self) -> P2D {
        P2D {
            x: self.x.abs(),
            y: self.y.abs(),
        }
    }

    pub fn max(&self, other: &P2D) -> P2D {
        P2D {
            x: max_f32(self.x, other.x),
            y: max_f32(self.y, other.y),
        }
    }

    pub fn min(&self, other: &P2D) -> P2D {
        P2D {
            x: min_f32(self.x, other.x),
            y: min_f32(self.y, other.y),
        }
    }

    pub fn unitize(&self) -> P2D {
        let m = self.mag();
        (1.0 / m) * self
    }

    pub fn normal(&self) -> P2D {
        P2D {
            x: -1.0 * self.y,
            y: self.x,
        }
    }

    pub fn scalar_mul_x(&self, s: f32) -> P2D {
        P2D {
            x: self.x * s,
            y: self.y,
        }
    }

    pub fn scalar_mul_y(&self, s: f32) -> P2D {
        P2D {
            x: self.x,
            y: self.y * s,
        }
    }

    pub fn powi(&self, x: i32) -> P2D {
        P2D {
            x: self.x.powi(x),
            y: self.y.powi(x),
        }
    }
}

impl Add for P2D {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        P2D {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Add for &P2D {
    type Output = P2D;

    fn add(self, rhs: Self) -> Self::Output {
        P2D {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Div for P2D {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        P2D {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
        }
    }
}

impl Mul<P2D> for f32 {
    type Output = P2D;

    fn mul(self, rhs: P2D) -> Self::Output {
        P2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<&P2D> for f32 {
    type Output = P2D;

    fn mul(self, rhs: &P2D) -> Self::Output {
        P2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<P2D> for &f32 {
    type Output = P2D;

    fn mul(self, rhs: P2D) -> Self::Output {
        P2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Add<P2D> for f32 {
    type Output = P2D;

    fn add(self, rhs: P2D) -> Self::Output {
        P2D {
            x: self + rhs.x,
            y: self + rhs.y,
        }
    }
}

impl Sub for P2D {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        P2D {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<'a, 'b> Sub<&'b P2D> for &'a P2D {
    type Output = P2D;

    fn sub(self, other: &'b P2D) -> P2D {
        P2D {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Display for P2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

/// Determine four-quadrant arctan of (y/x) constrained between `0.0` and `2*PI`.
pub fn arctan(x: f32, y: f32) -> Radians {
    Radians(modulo_f32(y.atan2(x), RAD_2PI.0))
}

/// Calculate `x` modulo `y`.
pub fn modulo_f32(x: f32, y: f32) -> f32 {
    (x % y) + y
}

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

/// Calculate the area of an "ideal" initial cell of radius R, if it has n vertices.
pub fn calc_init_cell_area(r: f32, n: u32) -> f32 {
    let poly_coords = (0..n).map(|vix| {
        let theta = (vix as f32)/(n as f32) * 2.0 * PI;
        P2D {
            x: r * theta.cos(),
            y: r * theta.sin(),
        }
    }).collect::<Vec<P2D>>();
    calc_poly_area(&poly_coords)
}

// /// Given three points `p0`, `p1`, `p2`, check if `p2` is left of the line through `p0` and `p1`.
// /// Greater than 0 if `p2` is left of, `0` if `p2` is on, and less than 0 if `p2` is right of.
// pub fn is_point_left_of_line(p0: &P2D, p1: &P2D, p2: &P2D) -> f32 {
//     (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)
// }

// pub struct Bbox {
//     xmin: f32,
//     ymin: f32,
//     xmax: f32,
//     ymax: f32,
// }
//
// impl Bbox {
//     pub fn point_in(&self, p: &P2D) -> bool {
//         (self.xmin < p.x) && (self.ymin < p.y) && (self.xmax > p.x) && (self.ymax > p.y)
//     }
// }
//
// pub fn min_f32s(xs: &[f32]) -> f32 {
//     let mut r = -1.0 * std::f32::INFINITY;
//     for &x in xs {
//         if (x - r) < -1.0 * std::f32::EPSILON {
//             r = x
//         }
//     }
//     r
// }

pub fn min_f32(x: f32, y: f32) -> f32 {
    if x - y < 0.0 {
        x
    } else {
        y
    }
}

pub fn max_f32(x: f32, y: f32) -> f32 {
    if x - y > 0.0 {
        x
    } else {
        y
    }
}
//
// pub fn max_f32s(xs: &[f32]) -> f32 {
//     let mut r = std::f32::MIN;
//     for &x in xs {
//         r = max_f32(r, x);
//     }
//     r
// }

// pub fn poly_bbox(xys: &[P2D]) -> Bbox {
//     let xs: Vec<f32> = xys.iter().map(|v| v.x).collect();
//     let ys: Vec<f32> = xys.iter().map(|v| v.y).collect();
//     Bbox {
//         xmin: min_f32s(&xs),
//         ymin: min_f32s(&ys),
//         xmax: max_f32s(&xs),
//         ymax: max_f32s(&xs),
//     }
// }

// pub fn point_in_poly(xys: &[P2D], p: &P2D) -> bool {
//     if poly_bbox(xys).point_in(p) {
//         let nvs = xys.len();
//         let mut wn: i32 = 0;
//
//         (0..nvs).for_each(|i| {
//             let p_start = xys[i];
//             let p_end = xys[circ_ix_plus(i, nvs)];
//             let is_left = is_point_left_of_line(&p_start, &p_end, p);
//
//             if p_start.y <= p.y && p.y < p_end.y {
//                 // upward crossing
//                 if is_left > 0.0 {
//                     wn += 1;
//                 }
//             } else if p_end.y < p.y && p.y <= p_start.y {
//                 // downward crossing
//                 if is_point_left_of_line(&p_start, &p_end, p) < 0.0 {
//                     wn -= 1;
//                 }
//             }
//         });
//
//         if wn == 0 {
//             false
//         } else {
//             true
//         }
//     } else {
//         false
//     }
// }

pub fn hill_function3(thresh: f32, x: f32) -> f32 {
    let x_cubed = x.powi(3);
    x_cubed / (thresh.powi(3) + x_cubed)
}

/// Assumes that x, max_x > 0.
pub fn capped_linear_function(x: f32, max_x: f32) -> f32 {
    if x > max_x {
        1.0
    } else {
        x / max_x
    }
}
