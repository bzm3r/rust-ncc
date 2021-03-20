// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.inner <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.inner> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::close_to_zero;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::fmt;
use std::fmt::{Display, Formatter};
use std::ops::{Add, Mul, Sub};

#[derive(Copy, Clone, Default, Debug, Serialize, Deserialize)]
pub struct Radians {
    inner: f64,
}

impl From<f64> for Radians {
    fn from(val: f64) -> Self {
        Radians {
            inner: val.rem_euclid(2.0 * PI),
        }
    }
}

pub static RAD_0: Lazy<Radians> = Lazy::new(|| Radians::from(0.0));
pub static RAD_PI: Lazy<Radians> = Lazy::new(|| Radians::from(PI));
pub static RAD_2PI: Lazy<Radians> =
    Lazy::new(|| Radians::from(2.0 * PI));
pub static RAD_EPS: Lazy<Radians> = Lazy::new(|| Radians::from(1e-3));

impl PartialEq for Radians {
    fn eq(&self, other: &Self) -> bool {
        (self.inner - other.inner).abs() < RAD_EPS.inner
    }
}

/// Addition modulo `2*PI`.
impl Add for Radians {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Radians::from(self.inner + other.inner)
    }
}

/// Subtraction modulo `2*PI`.
impl Sub for Radians {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Radians::from(self.inner - other.inner)
    }
}

/// Addition modulo `2*PI`.
impl Add<&Radians> for Radians {
    type Output = Self;

    fn add(self, other: &Self) -> Self {
        Radians::from(self.inner + other.inner)
    }
}

/// Subtraction modulo `2*PI`.
impl Sub<&Radians> for Radians {
    type Output = Self;

    fn sub(self, other: &Self) -> Self {
        Radians::from(self.inner - other.inner)
    }
}

impl Mul<f64> for Radians {
    type Output = Radians;

    fn mul(self, rhs: f64) -> Self::Output {
        Radians::from(self.inner * rhs)
    }
}

impl Mul<Radians> for f64 {
    type Output = Radians;

    fn mul(self, rhs: Radians) -> Self::Output {
        Radians::from(self * rhs.inner)
    }
}

impl Radians {
    /// Calculates "circular" between.
    pub fn between(&self, t0: Radians, t1: Radians) -> bool {
        if close_to_zero(t0.inner - t1.inner, 1e-3) {
            false
        } else if t0.inner < t1.inner {
            t0.inner < self.inner && self.inner < t1.inner
        } else {
            !(t1.inner < self.inner && self.inner < t0.inner)
        }
    }
}

impl Display for Radians {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}_PI", self.inner / PI)
    }
}
