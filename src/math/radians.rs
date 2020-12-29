// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::modulo_f32;
use std::f32::consts::PI;
use std::ops::{Add, Sub};

/// Value always between `[0.0, 2*PI]`.
#[derive(
    PartialOrd, PartialEq, Clone, Copy,
)]
pub struct Radians(f32);

const RAD_2PI: Radians =
    Radians(2.0 * PI);
const RAD_EPS: Radians = Radians(1e-12);

/// Addition modulo `2*PI`.
impl Add for Radians {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Radians(modulo_f32(
            self.0 + other.0,
            RAD_2PI.0,
        ))
    }
}

/// Subtraction modulo `2*PI`.
impl Sub for Radians {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Radians(modulo_f32(
            self.0 - other.0,
            RAD_2PI.0,
        ))
    }
}

impl Radians {
    /// Helper function for `Radians::between`.
    fn _between(
        &self,
        t0: Radians,
        t1: Radians,
    ) -> bool {
        t0 < *self && *self < t1
    }

    /// Calculates if `self` is between `t0` and `t1`, where `t1` is assumed to be counter-clockwise after `t0`.
    pub fn between(
        &self,
        t0: Radians,
        t1: Radians,
    ) -> bool {
        if t0 - t1 < RAD_EPS {
            false
        } else if t1 - *self < RAD_EPS
            || t0 - *self < RAD_EPS
        {
            true
        } else if t0 < t1 {
            self._between(t0, t1)
        } else {
            !self._between(t1, t0)
        }
    }
}

/// Determine four-quadrant arctan of (y/x) constrained between `0.0` and `2*PI`.
pub fn arctan(
    x: f32,
    y: f32,
) -> Radians {
    Radians(modulo_f32(
        y.atan2(x),
        RAD_2PI.0,
    ))
}
