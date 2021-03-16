// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::exp_setup::exp_parser::{ExperimentArgs, RgtpDistribDefs};
use crate::math::v2d::V2d;
use crate::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::world::IntegratorOpts;
use crate::Directories;

pub mod defaults;
pub mod exp_parser;
pub mod markers;
pub mod n_cells;
pub mod pair;
pub mod py_compare;

use crate::parameters::quantity::Time;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum ExperimentType {
    NCells {
        num_cells: usize,
    },
    Pair {
        sep_in_cell_diams: usize,
        rgtp_distrib_defs_per_cell: Vec<RgtpDistribDefs>,
    },
    PyCompare {
        num_cells: usize,
        py_main: Option<PathBuf>,
    },
}

/// Generate the experiment, so that it can be run.
#[allow(clippy::too_many_arguments)]
pub fn generate(
    dirs: Directories,
    args: ExperimentArgs,
) -> Vec<Experiment> {
    dirs.make();
    match &args.ty {
        ExperimentType::NCells { .. } => {
            n_cells::generate(dirs, args)
        }
        ExperimentType::Pair { .. } => pair::generate(dirs, args),
        ExperimentType::PyCompare { .. } => {
            py_compare::generate(dirs, args)
        }
    }
}

/// Specifies initial placement of the group.
#[derive(Clone)]
pub struct GroupBBox {
    /// Width of group in terms of cell diameter.
    pub width: usize,
    /// Height of group in terms of cell diameter.
    pub height: usize,
    /// Bottom left of the group in normalized space units.
    pub bottom_left: V2d,
}

/// Information required for a cell group to be created.
#[derive(Clone)]
pub struct CellGroup {
    /// The number of cells in the group.
    pub num_cells: usize,
    /// Initial layout of the cell group.
    pub layout: GroupBBox,
    /// Parameters shared by all cells in this group.
    pub parameters: Parameters,
}

/// Information required to create an experiment.
#[derive(Clone)]
pub struct Experiment {
    pub ty: ExperimentType,
    pub name: String,
    /// Time period in seconds to be simulated.
    pub final_t: Time,
    /// Characteristic quantities.
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    /// List of cell groups involved in this experiment.
    pub cell_groups: Vec<CellGroup>,
    /// Random number generator to be used for various purposes.
    /// Initialized from a seed, otherwise from "entropy".
    pub rng: Pcg32,
    /// Seed that was used to initialize rng, if it generated from a
    /// seed.
    pub seed: u64,
    pub snap_period: Time,
    pub max_on_ram: usize,
    pub int_opts: IntegratorOpts,
    pub out_dir: PathBuf,
    pub py_main: Option<PathBuf>,
}
