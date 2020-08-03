// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use serde::{Deserialize, Serialize};
use avro_schema_derive::Schematize;
use crate::math::radians::{Radians, arctan};
use crate::math::{max_f32, min_f32};
use std::ops::{Add, Div, Mul, Sub};
use std::fmt::Display;
use std::fmt;

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
