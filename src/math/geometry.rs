// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::p2d::P2D;
use crate::utils::{circ_ix_minus, circ_ix_plus};

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

// /// Given three points `p0`, `p1`, `p2`, check if `p2` is left of the line through `p0` and `p1`.
// /// Greater than 0 if `p2` is left of, `0` if `p2` is on, and less than 0 if `p2` is right of.
// pub fn is_point_left_of_line(p0: &P2D, p1: &P2D, p2: &P2D) -> f32 {
//     (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y)
// }

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
