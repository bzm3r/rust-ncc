// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod single;

use crate::math::v2d::V2d;
use crate::parameters::{CharQuantities, Parameters, WorldParameters};

/// Specifies initial placement of the group.
pub struct GroupLayout {
    /// Width of group in terms of number of cells.
    pub width: u32,
    /// Height of group in terms of number if cells.
    pub height: u32,
    /// Bottom left of the group in micrometers.
    pub bottom_left: V2d,
}

/// Information required for a cell group to be created.
pub struct CellGroup {
    /// The number of cells in the group.
    pub num_cells: u32,
    /// Initial layout of the cell group.
    pub layout: GroupLayout,
    /// Parameters shared by all cells in this group.
    pub parameters: Parameters,
}

/// Information required to create an experiment.
pub struct Experiment {
    pub title: String,
    /// Characteristic quantities.
    pub char_quants: CharQuantities,
    pub world_parameters: WorldParameters,
    /// List of cell groups involved in this experiment.
    pub cell_groups: Vec<CellGroup>,
}
