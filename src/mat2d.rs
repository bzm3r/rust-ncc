// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

// pub struct Mat2D {
//     dat: Vec<f32>,
//     nr: u32,
//     nc: u32,
// }
//
// impl Mat2D {
//     pub fn new(nr: u32, nc: u32) -> Mat2D {
//         Mat2D {
//             dat: vec![0.0_f32; nr * nc],
//             nr,
//             nc,
//         }
//     }
//
//     pub fn get_row(&self, ix: u32) -> &[f32] {
//         let start = ix * self.nc;
//         let end = start + self.nc;
//         &self.dat[start..end]
//     }
//
//     pub fn set_row(&mut self, ix: u32, vs: &[f32]) {
//         let start = ix * self.nc;
//         let end = start + self.nc;
//         for (d, v) in self.dat[start..end].iter_mut().zip(vs) {
//             *d = *v;
//         }
//     }
//
//     pub fn get_row_mut(&mut self, ix: u32) -> &mut [f32] {
//         let start = ix * self.nc;
//         let end = start + self.nc;
//         &mut self.dat[start..end]
//     }
//
//     pub fn get_col(&self, ix: u32) -> Vec<f32> {
//         (0..self.nr)
//             .map(|i| self.dat[i * self.nr + ix])
//             .collect()
//     }
//
//     pub fn get(&self, x: u32, y: u32) -> f32 {
//         self.dat[x * self.nc + y]
//     }
//
//     pub fn set(&mut self, x: u32, y: u32, v: f32) {
//         self.dat[x * self.nc + y] = v;
//     }
// }
