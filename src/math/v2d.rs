// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::radians::{arctan, Radians};
use crate::math::{max_f32, min_f32};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::fmt::Write;
use std::ops::{Add, Div, Mul, Sub};

/// 2D vector with `f32` elements.
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct V2D {
    pub x: f32,
    pub y: f32,
}

impl V2D {
    /// Create a new `V2d` given two coordinate points `(x, y)` of
    /// type `(f32, f32)`.
    pub fn new(x: f32, y: f32) -> V2D {
        V2D { x, y }
    }

    /// Calculate direction of vector from origin to point in radians.
    pub fn direction(&self) -> Radians {
        arctan(self.x, self.y)
    }

    /// Calculate magnitude of vector from origin to point.
    #[inline]
    pub fn mag(&self) -> f32 {
        self.mag_squared().sqrt()
    }

    /// Calculate magnitude squared of vector.
    #[inline]
    pub fn mag_squared(&self) -> f32 {
        self.x.powi(2) + self.y.powi(2)
    }

    pub fn dot(&self, other: &V2D) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn abs(&self) -> V2D {
        V2D {
            x: self.x.abs(),
            y: self.y.abs(),
        }
    }

    pub fn max(&self, other: &V2D) -> V2D {
        V2D {
            x: max_f32(self.x, other.x),
            y: max_f32(self.y, other.y),
        }
    }

    pub fn min(&self, other: &V2D) -> V2D {
        V2D {
            x: min_f32(self.x, other.x),
            y: min_f32(self.y, other.y),
        }
    }

    pub fn unitize(&self) -> V2D {
        let m = self.mag();
        (1.0 / m) * self
    }

    pub fn normal(&self) -> V2D {
        V2D {
            x: -1.0 * self.y,
            y: self.x,
        }
    }

    pub fn scalar_mul_x(&self, s: f32) -> V2D {
        V2D {
            x: self.x * s,
            y: self.y,
        }
    }

    pub fn scalar_mul_y(&self, s: f32) -> V2D {
        V2D {
            x: self.x,
            y: self.y * s,
        }
    }

    pub fn powi(&self, x: i32) -> V2D {
        V2D {
            x: self.x.powi(x),
            y: self.y.powi(x),
        }
    }

    pub fn scale(&self, factor: f32) -> V2D {
        V2D {
            x: self.x * factor,
            y: self.y * factor,
        }
    }

    pub fn translate(&self, dx: f32, dy: f32) -> V2D {
        V2D {
            x: self.x + dx,
            y: self.y + dy,
        }
    }

    pub fn close_to_zero(&self) -> bool {
        self.x.abs() < f32::EPSILON && self.y.abs() < f32::EPSILON
    }

    #[inline]
    pub fn zeros() -> V2D {
        V2D::default()
    }
}

impl Add for V2D {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        V2D {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Add for &V2D {
    type Output = V2D;

    fn add(self, rhs: Self) -> Self::Output {
        V2D {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Div for V2D {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        V2D {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
        }
    }
}

impl Mul<V2D> for f32 {
    type Output = V2D;

    fn mul(self, rhs: V2D) -> Self::Output {
        V2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<&V2D> for f32 {
    type Output = V2D;

    fn mul(self, rhs: &V2D) -> Self::Output {
        V2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Mul<V2D> for &f32 {
    type Output = V2D;

    fn mul(self, rhs: V2D) -> Self::Output {
        V2D {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Add<V2D> for f32 {
    type Output = V2D;

    fn add(self, rhs: V2D) -> Self::Output {
        V2D {
            x: self + rhs.x,
            y: self + rhs.y,
        }
    }
}

impl Sub for V2D {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        V2D {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<'a, 'b> Sub<&'b V2D> for &'a V2D {
    type Output = V2D;

    fn sub(self, other: &'b V2D) -> V2D {
        V2D {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Display for V2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

#[allow(unused)]
pub fn poly_to_string(poly: &[V2D]) -> String {
    let mut r = "[".to_string();
    for p in poly {
        writeln!(r, "{},", p).unwrap();
    }
    write!(r, "]").unwrap();
    r
}
