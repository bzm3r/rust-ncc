// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub fn circ_ix_minus(ix: usize, num_items: usize) -> usize {
    if ix == 0 {
        num_items - 1
    } else {
        ix - 1
    }
}

pub fn circ_ix_plus(ix: usize, num_items: usize) -> usize {
    if ix == (num_items - 1) {
        0
    } else {
        ix + 1
    }
}

pub fn circ_ix(ix: usize, delta: isize, num_items: usize) -> usize {
    let q = delta / (num_items as isize);
    let delta = delta - q * (num_items as isize);

    let shift_x = (ix as isize) + delta;
    if shift_x < 0 {
        num_items - (shift_x.abs() as usize)
    } else {
        (shift_x as usize) - num_items
    }
}
