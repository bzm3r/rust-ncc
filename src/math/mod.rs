// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod geometry;
pub mod p2d;
pub mod radians;

/// Calculate `x` modulo `y`.
pub fn modulo_f32(x: f32, y: f32) -> f32 {
    (x % y) + y
}

/// Returns the minimum of a slice of floats.
pub fn min_f32s(xs: &[f32]) -> f32 {
    let mut x = std::f32::MAX;
    for &y in xs {
        if y < x {
            x = y;
        }
    }
    x
}

/// Returns the maximum of a slice of floats.
pub fn max_f32s(xs: &[f32]) -> f32 {
    let mut x = std::f32::MIN;
    for &y in xs {
        if y > x {
            x = y;
        }
    }
    x
}

/// Returns the minimum of two floats.
pub fn min_f32(x: f32, y: f32) -> f32 {
    if x < y {
        x
    } else {
        y
    }
}

/// Returns maximum of two floats.
pub fn max_f32(x: f32, y: f32) -> f32 {
    if x > y {
        x
    } else {
        y
    }
}

/// Hill function with exponent 3.
pub fn hill_function3(thresh: f32, x: f32) -> f32 {
    let x_cubed = x.powi(3);
    x_cubed / (thresh.powi(3) + x_cubed)
}

/// If `x > max_x`, returns 1.0, else `x/max_x`.
pub fn capped_linear_function(x: f32, max_x: f32) -> f32 {
    if x > max_x {
        1.0
    } else {
        x / max_x
    }
}
