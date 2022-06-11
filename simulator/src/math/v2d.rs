// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::radians::Radians;
use crate::math::{close_to_zero, max_f64, min_f64};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::fmt::Write;
use std::ops::{Add, Div, Mul, Sub};

/// 2D vector with `f64` elements.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq)]
pub struct V2d {
    pub x: f64,
    pub y: f64,
}

impl V2d {
    /// Create a new `V2d` given two coordinate points `(x, y)` of
    /// type `(f64, f64)`.
    pub fn new(x: f64, y: f64) -> V2d {
        V2d { x, y }
    }

    /// Calculate direction of vector from origin to point in radians.
    pub fn direction(&self) -> Radians {
        Radians::from(self.y.atan2(self.x))
    }

    /// Calculate magnitude of vector from origin to point.
    #[inline]
    pub fn mag(&self) -> f64 {
        self.mag_squared().sqrt()
    }

    /// Calculate magnitude squared of vector.
    #[inline]
    pub fn mag_squared(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }

    pub fn dot(&self, other: &V2d) -> f64 {
        self.x * other.x + self.y * other.y
    }

    pub fn abs(&self) -> V2d {
        V2d {
            x: self.x.abs(),
            y: self.y.abs(),
        }
    }

    pub fn max(&self, other: &V2d) -> V2d {
        V2d {
            x: max_f64(self.x, other.x),
            y: max_f64(self.y, other.y),
        }
    }

    pub fn min(&self, other: &V2d) -> V2d {
        V2d {
            x: min_f64(self.x, other.x),
            y: min_f64(self.y, other.y),
        }
    }

    pub fn unitize(&self) -> V2d {
        let m = self.mag();
        (1.0 / m) * self
    }

    pub fn normal(&self) -> V2d {
        V2d {
            x: -1.0 * self.y,
            y: self.x,
        }
    }

    pub fn scalar_mul_x(&self, s: f64) -> V2d {
        V2d {
            x: self.x * s,
            y: self.y,
        }
    }

    pub fn scalar_mul_y(&self, s: f64) -> V2d {
        V2d {
            x: self.x,
            y: self.y * s,
        }
    }

    pub fn powi(&self, x: i32) -> V2d {
        V2d {
            x: self.x.powi(x),
            y: self.y.powi(x),
        }
    }

    pub fn scale(&self, factor: f64) -> V2d {
        V2d {
            x: self.x * factor,
            y: self.y * factor,
        }
    }

    pub fn translate(&self, dx: f64, dy: f64) -> V2d {
        V2d {
            x: self.x + dx,
            y: self.y + dy,
        }
    }

    pub fn close_to_zero(&self) -> bool {
        self.x.abs() < f64::EPSILON && self.y.abs() < f64::EPSILON
    }

    #[inline]
    pub fn zeros() -> V2d {
        V2d::default()
    }

    pub fn close_to(&self, other: &V2d, eps: f64) -> bool {
        close_to_zero(self.x - other.x, eps)
            && close_to_zero(self.y - other.y, eps)
    }

    pub fn cross(&self, other: &V2d) -> f64 {
        self.x * other.y - self.y * other.x
    }
}

impl Add for V2d {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        V2d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Add for &V2d {
    type Output = V2d;

    fn add(self, rhs: Self) -> Self::Output {
        V2d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Add<f64> for V2d {
    type Output = Self;

    fn add(self, rhs: f64) -> Self::Output {
        V2d {
            x: self.x + rhs,
            y: self.y + rhs,
        }
    }
}

impl Add<V2d> for f64 {
    type Output = V2d;

    fn add(self, rhs: V2d) -> Self::Output {
        rhs + self
    }
}

impl Div for V2d {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        V2d {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
        }
    }
}

impl Mul<V2d> for f64 {
    type Output = V2d;

    fn mul(self, rhs: V2d) -> Self::Output {
        V2d {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<&V2d> for f64 {
    type Output = V2d;

    fn mul(self, rhs: &V2d) -> Self::Output {
        V2d {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<V2d> for &f64 {
    type Output = V2d;

    fn mul(self, rhs: V2d) -> Self::Output {
        V2d {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Sub for V2d {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        V2d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<'a, 'b> Sub<&'b V2d> for &'a V2d {
    type Output = V2d;

    fn sub(self, other: &'b V2d) -> V2d {
        V2d {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Display for V2d {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

pub fn poly_to_string(poly: &[V2d]) -> String {
    let mut r = "[".to_string();
    for p in poly {
        writeln!(r, "{},", p).unwrap();
    }
    write!(r, "]").unwrap();
    r
}

/// Square of a 2D vector with `f64` elements.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq)]
pub struct SqP2d {
    pub x: f64,
    pub y: f64,
}

impl Mul<V2d> for V2d {
    type Output = SqP2d;

    fn mul(self, rhs: Self) -> Self::Output {
        SqP2d {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
        }
    }
}

impl Div<V2d> for SqP2d {
    type Output = V2d;

    fn div(self, rhs: V2d) -> Self::Output {
        V2d {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
        }
    }
}
