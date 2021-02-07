// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod geometry;
pub mod radians;
pub mod v2d;

#[inline]
/// Calculate `x` modulo `y`.
pub fn modulo_f64(x: f64, y: f64) -> f64 {
    (x % y) + y
}

/// Returns the minimum of a slice of floats.
pub fn min_f64s(xs: &[f64]) -> f64 {
    let mut x = std::f64::MAX;
    for &y in xs {
        if y < x {
            x = y;
        }
    }
    x
}

/// Returns the maximum of a slice of floats.
pub fn max_f64s(xs: &[f64]) -> f64 {
    let mut x = std::f64::MIN;
    for &y in xs {
        if y > x {
            x = y;
        }
    }
    x
}

#[inline]
/// Returns the minimum of two floats.
pub fn min_f64(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

#[inline]
/// Returns maximum of two floats.
pub fn max_f64(x: f64, y: f64) -> f64 {
    if x > y {
        x
    } else {
        y
    }
}

#[inline]
/// Hill function with exponent 3.
pub fn hill_function3(thresh: f64, x: f64) -> f64 {
    let x_cubed = x.powi(3);
    x_cubed / (thresh.powi(3) + x_cubed)
}

#[inline]
/// Model a function which is `0` and `1` at the given inputs,
/// but is capped in output between `[0, 1]`.
pub fn capped_linear_fn(x: f64, zero_at: f64, one_at: f64) -> f64 {
    // first, calculate the uncapped linear function
    let m = 1.0 / (one_at - zero_at);
    // 0.0 = y = m * zero_at + b => 0.0 - m * zero_at = b
    let b = -1.0 * m * zero_at;
    let y = m * x + b;
    if y < 0.0 {
        0.0
    } else if y > 1.0 {
        1.0
    } else {
        y
    }
}

/// Return if the float `x` close to `0.0`.
pub fn close_to_zero(x: f64) -> bool {
    x.abs() < 1e-4
}

/// Round to `n` digits.
pub fn round(x: f64, n: u32) -> f64 {
    let scale = 10_u32.pow(n) as f64;
    (x * scale).round() / scale
}

#[derive(Clone, Copy)]
pub enum InUnitInterval {
    In,
    Zero,
    One,
    Out,
}

#[inline]
pub fn in_unit_interval(x: f64) -> InUnitInterval {
    if x > 0.0 && x < 1.0 {
        InUnitInterval::In
    } else if close_to_zero(x) {
        InUnitInterval::Zero
    } else if close_to_zero(1.0 - x) {
        InUnitInterval::One
    } else {
        InUnitInterval::Out
    }
}
