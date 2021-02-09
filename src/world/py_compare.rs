// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::cell::py_compare::{Cell, PrintOptions};
use crate::experiments::Experiment;
use crate::hardio::py_compare::Writer;
use crate::interactions::{
    InteractionGenerator, Interactions, RelativeRgtpActivity,
};
use crate::math::v2d::V2D;
use crate::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::NVERTS;

use crate::cell::states::Core;
use crate::world::gen_cell_centroids;
use rand::{RngCore, SeedableRng};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Clone, Deserialize, Serialize, PartialEq, Default, Debug)]
pub struct Cells {
    pub states: Vec<Cell>,
    pub interactions: Vec<Interactions>,
}

impl Cells {
    #[allow(clippy::too_many_arguments)]
    fn simulate(
        &self,
        tstep: u32,
        num_int_steps: u32,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
        writer: &mut Writer,
    ) -> Result<Cells, String> {
        let mut new_cell_states =
            vec![self.states[0]; self.states.len()];
        let shuffled_cells = {
            let crs = self.states.iter().collect::<Vec<&Cell>>();
            // crs.shuffle(rng);
            crs
        };
        for cell_state in shuffled_cells {
            let ci = cell_state.ix;
            let contact_data =
                interaction_generator.get_contact_data(ci);

            let new_cell_state = cell_state.simulate_euler(
                tstep,
                num_int_steps,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cell_state.group_ix],
                rng,
                writer,
            )?;
            interaction_generator
                .update(ci, &new_cell_state.core.poly);

            new_cell_states[ci] = new_cell_state;
        }
        let rel_rgtps = new_cell_states
            .iter()
            .map(|c| {
                c.core.calc_relative_rgtp_activity(
                    &group_parameters[c.group_ix],
                )
            })
            .collect::<Vec<[RelativeRgtpActivity; NVERTS]>>();
        Ok(Cells {
            states: new_cell_states,
            interactions: interaction_generator.generate(&rel_rgtps),
        })
    }
}

#[derive(Clone, Deserialize, Serialize, PartialEq, Debug)]
pub struct Snapshot {
    pub tstep: u32,
    pub cells: Cells,
    pub rng: Pcg32,
}

#[derive(Deserialize, Serialize, Clone, Default, Debug, PartialEq)]
pub struct WorldInfo {
    pub snap_freq: u32,
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    pub cell_params: Vec<Parameters>,
}

impl Iterator for Cells {
    type Item = ();

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

pub struct World {
    char_quants: CharQuantities,
    tstep: u32,
    world_params: WorldParameters,
    group_params: Vec<Parameters>,
    cells: Cells,
    interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    writer: Writer,
}

fn gen_poly(centroid: &V2D, radius: f64) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    (0..NVERTS).for_each(|vix| {
        let vf = (vix as f64) / (NVERTS as f64);
        let theta = 2.0 * PI * vf;
        r[vix] = V2D {
            x: centroid.x + theta.cos() * radius,
            y: centroid.y + theta.sin() * radius,
        };
    });
    r
}

impl World {
    pub fn new(experiment: Experiment) -> World {
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
            .map(|cg| cg.parameters)
            .collect::<Vec<Parameters>>();

        // Create a list of indices of the groups and create a vector
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
                Core::init(
                    *poly,
                    parameters.init_rac,
                    parameters.init_rho,
                )
            })
            .collect::<Vec<Core>>();
        // Calculate relative activity of Rac1 vs. RhoA at a node.
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
        let cell_interactions =
            interaction_generator.generate(&cell_rgtps);
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
                PrintOptions {
                    deltas: false,
                    interaction_updates: false,
                },
            ));
        }
        let cells = Cells {
            states: cell_states,
            interactions: cell_interactions,
        };
        World {
            tstep: 0,
            char_quants,
            world_params,
            group_params,
            cells,
            interaction_generator,
            rng,
            writer: Writer::default(),
        }
    }

    pub fn take_snapshot(&self) -> Snapshot {
        Snapshot {
            tstep: self.tstep,
            rng: self.rng,
            cells: self.cells.clone(),
        }
    }

    pub fn simulate(&mut self, final_tpoint: f64) {
        let num_tsteps =
            (final_tpoint / self.char_quants.time()).ceil() as u32;
        let num_int_steps = 10;
        let mut writer = self.writer.init(
            num_tsteps as usize,
            num_int_steps as usize,
            self.cells.states.len(),
            self.world_params.interactions.phys_contact.cil_mag
                as u32,
            self.world_params.interactions.coa.map_or_else(
                || 0,
                |p| (p.vertex_mag as f64 * NVERTS as f64) as u32,
            ),
        );
        writer.save_header(
            &self.char_quants,
            &self.world_params,
            &self.group_params[0],
        );
        while self.tstep < num_tsteps {
            let new_cells: Cells = self
                .cells
                .simulate(
                    self.tstep,
                    num_int_steps,
                    &mut self.rng,
                    &self.world_params,
                    &self.group_params,
                    &mut self.interaction_generator,
                    &mut writer,
                )
                .unwrap();

            self.cells = new_cells;
            self.tstep += 1;
        }
        if !writer.finished {
            panic!(
                format!(
                    "Unfinished writer. num_tsteps to store: {}, num_tsteps stored: {}",
                    writer.num_tsteps, writer.data.tsteps.len()
                )
            )
        }
    }
}
