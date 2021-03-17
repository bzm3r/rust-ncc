// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::exp_setup::exp_parser::ExperimentArgs;
use crate::math::v2d::V2d;
use crate::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::world::IntegratorOpts;
use crate::{Directories, NVERTS};

pub mod defaults;
pub mod exp_parser;
pub mod markers;
pub mod n_cells;
pub mod pair;
pub mod py_compare;

use crate::cell::chemistry::distrib_gens::{
    random, specific_random, specific_uniform,
};
use crate::exp_setup::markers::mark_verts;
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
        rgtp_distrib_defs_per_cell: PairRgtpDistribDefs,
    },
    PyCompare {
        num_cells: usize,
        py_main: Option<PathBuf>,
    },
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum DistribDef {
    Random { frac: f64 },
    SpecificRandom { frac: f64, marked_verts: Vec<usize> },
    SpecificUniform { frac: f64, marked_verts: Vec<usize> },
}

impl DistribDef {
    pub fn into_distrib(&self, rng: &mut Pcg32) -> [f64; NVERTS] {
        match self {
            DistribDef::Random { frac } => random(rng, *frac),
            DistribDef::SpecificRandom { frac, marked_verts } => {
                specific_random(rng, *frac, mark_verts(&marked_verts))
            }
            DistribDef::SpecificUniform { frac, marked_verts } => {
                specific_uniform(*frac, mark_verts(&marked_verts))
            }
        }
    }
}

impl Default for DistribDef {
    fn default() -> Self {
        DistribDef::Random { frac: 0.1 }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RgtpDistribDef {
    acts: DistribDef,
    inacts: DistribDef,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct RgtpDistribDefs {
    rac: RgtpDistribDef,
    rho: RgtpDistribDef,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
pub struct PairRgtpDistribDefs {
    cell0: RgtpDistribDefs,
    cell1: RgtpDistribDefs,
}

impl Default for ExperimentType {
    fn default() -> Self {
        ExperimentType::NCells { num_cells: 1 }
    }
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
