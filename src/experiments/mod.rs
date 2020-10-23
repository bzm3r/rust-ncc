// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod adh;
pub mod adh2;
pub mod cil;
pub mod single;

use crate::math::p2d::V2D;
use crate::parameters::{BasicQuants, GlobalParameters, Parameters};

pub struct GroupLayout {
    pub width: u32,
    pub height: u32,
    pub bottom_left: V2D,
}

pub struct CellGroup {
    pub num_cells: u32,
    pub layout: GroupLayout,
    pub parameters: Parameters,
}

pub struct Experiment {
    pub title: String,
    pub basic_quants: BasicQuants,
    pub world_parameters: GlobalParameters,
    pub cell_groups: Vec<CellGroup>,
}
