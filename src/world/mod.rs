// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod hardio;
#[cfg(feature = "validate")]
use crate::cell::confirm_volume_exclusion;
use crate::cell::core_state::CoreState;
use crate::cell::Cell;
use crate::experiments::{CellGroup, Experiment};
use crate::interactions::{
    InteractionGenerator, Interactions, RelativeRgtpActivity,
};
use crate::math::v2d::V2D;
use crate::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use crate::NVERTS;
//use rand_core::SeedableRng;
use crate::utils::pcg32::Pcg32;
use crate::world::hardio::{save_compact, save_full, Format};
use rand::seq::SliceRandom;
use rand::{RngCore, SeedableRng};
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::f32::consts::PI;
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize)]
pub struct Cells {
    pub states: Vec<Cell>,
    pub interactions: Vec<Interactions>,
}

impl Cells {
    fn simulate(
        &self,
        tstep: u32,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
    ) -> Result<Cells, String> {
        let mut new_cell_states =
            vec![self.states[0]; self.states.len()];
        let shuffled_cells = {
            let mut crs = self.states.iter().collect::<Vec<&Cell>>();
            crs.shuffle(rng);
            crs
        };
        for cell_state in shuffled_cells {
            let ci = cell_state.ix;
            let contact_data =
                interaction_generator.get_contact_data(ci);

            let new_cell_state = cell_state.simulate_rkdp5(
                tstep,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cell_state.group_ix],
                rng,
            )?;

            interaction_generator
                .update(ci, &new_cell_state.core.poly);

            #[cfg(feature = "validate")]
            confirm_volume_exclusion(
                &new_cell_state.core.poly,
                &interaction_generator.get_contact_data(ci),
                &format!("world cell {}", ci),
            )?;

            new_cell_states[ci] = new_cell_state;
        }
        Ok(Cells {
            states: new_cell_states,
            interactions: interaction_generator.generate(),
        })
    }
}

#[derive(Deserialize, Serialize, Clone)]
pub struct DeepSnapshot {
    pub tstep: u32,
    pub interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    pub cells: Cells,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct Snapshot {
    pub tstep: u32,
    pub cells: Cells,
}

impl DeepSnapshot {
    pub fn to_mini(&self) -> Snapshot {
        Snapshot {
            tstep: self.tstep,
            cells: self.cells.clone(),
        }
    }
}

#[derive(Deserialize, Serialize, Clone)]
pub struct History {
    pub snap_freq: u32,
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    pub cell_params: Vec<Parameters>,
    pub snapshots: Vec<Snapshot>,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct DeepHistory {
    pub snap_freq: u32,
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    pub cell_params: Vec<Parameters>,
    pub snapshots: Vec<DeepSnapshot>,
}

impl DeepHistory {
    pub fn to_history(&self) -> History {
        History {
            snap_freq: self.snap_freq,
            char_quants: self.char_quants,
            world_params: self.world_params.clone(),
            cell_params: self.cell_params.clone(),
            snapshots: self
                .snapshots
                .iter()
                .map(|fs| fs.to_mini())
                .collect::<Vec<Snapshot>>(),
        }
    }
}

pub struct World {
    char_quants: CharQuantities,
    tstep: u32,
    world_params: WorldParameters,
    group_params: Vec<Parameters>,
    pub history: Vec<DeepSnapshot>,
    cells: Cells,
    interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    out_dir: PathBuf,
    file_name: String,
    snap_freq: u32,
}

fn gen_poly(centroid: &V2D, radius: f32) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    (0..NVERTS).for_each(|vix| {
        let vf = (vix as f32) / (NVERTS as f32);
        let theta = 2.0 * PI * vf;
        r[vix] = V2D {
            x: centroid.x + theta.cos() * radius,
            y: centroid.y + theta.sin() * radius,
        };
    });
    r
}

impl World {
    pub fn new(
        experiment: Experiment,
        output_dir: PathBuf,
        snap_freq: u32,
    ) -> World {
        // Unpack relevant info from `Experiment` data structure.
        let Experiment {
            char_quants,
            world_parameters: world_params,
            cell_groups,
            mut rng,
            ..
        } = experiment;
        // Extract the parameters from each `CellGroup` object obtained
        // from the `Experiment`.
        let group_params = cell_groups
            .iter()
            .map(|cg| cg.parameters.clone())
            .collect::<Vec<Parameters>>();

        // Create a list of indices of the groups. and reate a vector
        // of the cell centroids in each group.
        let mut cell_group_ixs = vec![];
        let mut cell_centroids = vec![];
        cell_groups.iter().enumerate().for_each(|(gix, cg)| {
            cell_group_ixs.append(&mut vec![gix; cg.num_cells]);
            cell_centroids
                .append(&mut gen_cell_centroids(cg).unwrap())
        });
        // Generate the cell polygons from the cell centroid
        // information generated in the last step.
        let cell_polys = cell_group_ixs
            .iter()
            .zip(cell_centroids.iter())
            .map(|(&gix, cc)| gen_poly(cc, group_params[gix].cell_r))
            .collect::<Vec<[V2D; NVERTS]>>();
        // Create initial cell states, using the parameters associated
        // with the cell a cell group is in, and the cell's centroid
        // location.
        let cell_core_states = cell_group_ixs
            .iter()
            .zip(cell_polys.iter())
            .map(|(&gix, poly)| {
                let parameters = &group_params[gix];
                CoreState::new(
                    *poly,
                    parameters.init_rac,
                    parameters.init_rho,
                )
            })
            .collect::<Vec<CoreState>>();
        // Calcualte relative activity of Rac1 vs. RhoA at a node.
        // This is needed for CRL.
        let cell_rgtps = cell_group_ixs
            .iter()
            .zip(cell_core_states.iter())
            .map(|(&gix, state)| {
                let parameters = &group_params[gix];
                state.calc_relative_rgtp_activity(parameters)
            })
            .collect::<Vec<[RelativeRgtpActivity; NVERTS]>>();
        // Create a new `InteractionGenerator`.
        let interaction_generator = InteractionGenerator::new(
            &cell_polys,
            &cell_rgtps,
            world_params.interactions.clone(),
        );
        // Generate initial cell interactions.
        let cell_interactions = interaction_generator.generate();
        // Create `Cell` structures to represent each cell, and the random number generator associated per cell.
        let mut cell_states = vec![];
        for (cell_ix, group_ix) in
            cell_group_ixs.into_iter().enumerate()
        {
            // Parameters that will be used by this cell. Determined
            // by figuring out which group it belongs to, as all cells
            // within a group use the same parameters.
            let parameters = &group_params[group_ix];
            let mut cell_rng = Pcg32::seed_from_u64(rng.next_u64());
            // Create a new cell.
            cell_states.push(Cell::new(
                cell_ix,
                group_ix,
                cell_core_states[cell_ix],
                &cell_interactions[cell_ix],
                parameters,
                &mut cell_rng,
            ));
        }
        let cells = Cells {
            states: cell_states,
            interactions: cell_interactions.clone(),
        };
        let history = vec![DeepSnapshot {
            tstep: 0,
            interaction_generator: interaction_generator.clone(),
            rng,
            cells: cells.clone(),
        }];
        World {
            tstep: 0,
            char_quants,
            world_params,
            group_params,
            history,
            cells,
            interaction_generator,
            rng,
            out_dir: output_dir,
            file_name: experiment.file_name,
            snap_freq,
        }
    }

    pub fn take_snapshot(&self) -> DeepSnapshot {
        DeepSnapshot {
            tstep: self.tstep,
            interaction_generator: self.interaction_generator.clone(),
            rng: self.rng,
            cells: self.cells.clone(),
        }
    }

    pub fn simulate(&mut self, final_tpoint: f32) {
        let num_tsteps =
            (final_tpoint / self.char_quants.time()).ceil() as u32;
        while self.tstep < num_tsteps {
            let new_cells: Cells = self
                .cells
                .simulate(
                    self.tstep,
                    &mut self.rng,
                    &self.world_params,
                    &self.group_params,
                    &mut self.interaction_generator,
                )
                .unwrap_or_else(|e| {
                    self.save_history(
                        true,
                        vec![Format::Cbor, Format::Bincode],
                    )
                    .unwrap();
                    panic!("tstep: {}\n{}", self.tstep, e);
                });

            self.cells = new_cells;
            self.tstep += 1;
            if self.tstep % self.snap_freq == 0 {
                self.history.push(self.take_snapshot());
            }
        }
    }

    pub fn deep_history(&self) -> DeepHistory {
        DeepHistory {
            snap_freq: self.snap_freq,
            char_quants: self.char_quants,
            world_params: self.world_params.clone(),
            cell_params: self
                .cells
                .states
                .iter()
                .map(|s| self.group_params[s.group_ix].clone())
                .collect::<Vec<Parameters>>(),
            snapshots: self.history.clone(),
        }
    }

    pub fn history(&self) -> History {
        History {
            snap_freq: self.snap_freq,
            char_quants: self.char_quants,
            world_params: self.world_params.clone(),
            cell_params: self
                .cells
                .states
                .iter()
                .map(|s| self.group_params[s.group_ix].clone())
                .collect::<Vec<Parameters>>(),
            snapshots: self
                .history
                .iter()
                .map(|fs| fs.to_mini())
                .collect::<Vec<Snapshot>>(),
        }
    }

    pub fn save_history(
        &self,
        compact: bool,
        formats: Vec<Format>,
    ) -> Result<(), Box<dyn Error>> {
        let history = self.history();
        println!("num snapshots: {}", history.snapshots.len());
        if compact {
            save_compact(
                history,
                &self.out_dir,
                formats,
                &self.file_name,
            )?;
        } else {
            save_full(
                self.deep_history(),
                &self.out_dir,
                formats,
                &self.file_name,
            )?;
        }
        Ok(())
    }
}

fn gen_cell_centroids(cg: &CellGroup) -> Result<Vec<V2D>, String> {
    let CellGroup {
        num_cells,
        layout,
        parameters,
    } = cg;
    let cell_r = parameters.cell_r;
    if layout.width * layout.height >= *num_cells {
        let mut r = vec![];
        let first_cell_centroid = V2D {
            x: layout.bottom_left.x + cell_r,
            y: layout.bottom_left.y + cell_r,
        };
        let row_delta = V2D {
            x: 0.0,
            y: 2.0 * cell_r,
        };
        let col_delta = V2D {
            x: 2.0 * cell_r,
            y: 0.0,
        };
        for ix in 0..*num_cells {
            let row = ix / layout.width;
            let col = ix - layout.width * row;
            let cg = first_cell_centroid
                + (row as f32) * row_delta
                + (col as f32) * col_delta;
            r.push(cg);
        }
        Ok(r)
    } else {
        Err(format!(
            "Cell group layout area ({}x{}) not large enough to fit {} cells.",
            layout.width, layout.height, num_cells
        ))
    }
}
