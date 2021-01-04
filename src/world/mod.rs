// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod hardio;
#[cfg(feature = "debug_mode")]
use crate::cell::confirm_volume_exclusion;
use crate::cell::core_state::CoreState;
use crate::cell::CellState;
use crate::experiments::{CellGroup, Experiment};
use crate::interactions::{
    CellInteractions, InteractionGenerator, RgtpActivityDiff,
};
use crate::math::v2d::V2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::NVERTS;
//use rand_core::SeedableRng;
use crate::utils::pcg32::Pcg32;
use rand::seq::SliceRandom;
use rand::{RngCore, SeedableRng};
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;
use std::fs::OpenOptions;
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize)]
pub struct Cells {
    pub cell_states: Vec<CellState>,
    pub interactions: Vec<CellInteractions>,
}

impl Cells {
    #[cfg(not(feature = "debug_mode"))]
    fn simulate(
        &self,
        rng: &mut Pcg32,
        cell_regs: &mut [Option<RandomEventGenerator>],
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
    ) -> Cells {
        let mut cells =
            vec![self.cell_states[0]; self.cell_states.len()];
        let shuffled_cells = {
            let mut crs =
                self.cell_states.iter().collect::<Vec<&CellState>>();
            crs.shuffle(rng);
            crs
        };
        for c in shuffled_cells {
            // println!("------------------");
            let ci = c.ix as usize;
            // println!("ci: {}", ci);
            // println!(
            //     "contacts: {:?}",
            //     interaction_generator.get_physical_contacts(ci)
            // );
            let contact_polys =
                interaction_generator.get_contact_data(ci);
            let new_cell_state = c.simulate_euler(
                self.tstep,
                &self.interactions[ci],
                contact_polys,
                cell_regs[ci].as_mut(),
                world_parameters,
                &group_parameters[c.group_ix as usize],
            );
            // println!("----------------------");
            interaction_generator
                .update(ci, &new_cell_state.state.vertex_coords);
            cells[ci] = new_cell_state;
        }
        Cells {
            tstep: self.tstep + 1,
            cell_states: cells,
            interactions: interaction_generator.generate(),
        }
    }

    #[cfg(feature = "debug_mode")]
    fn simulate(
        &self,
        tstep: u32,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
    ) -> Result<Cells, String> {
        dbg!("simulating cells...");
        let mut new_cell_states =
            vec![self.cell_states[0]; self.cell_states.len()];
        let shuffled_cells = {
            let mut crs =
                self.cell_states.iter().collect::<Vec<&CellState>>();
            crs.shuffle(rng);
            crs
        };
        for cell_state in shuffled_cells {
            // println!("------------------");
            let ci = cell_state.ix as usize;
            // println!("ci: {}", ci);
            // println!(
            //     "contacts: {:?}",
            //     interaction_generator.get_physical_contacts(ci)
            // );
            // let contacts =
            //     interaction_generator.get_physical_contacts(ci);
            let contact_data =
                interaction_generator.get_contact_data(ci);
            // println!("----------------------------");
            // println!(
            //     "    ci = {}; vi = {}; x_adh = {}; f_tot: {}",
            //     ci,
            //     vi,
            //     self.interactions[ci].x_adhs[vi].mag(),
            //     cell_state.mech.sum_fs[vi].mag(),
            // );
            let new_cell_state = cell_state.simulate_rkdp5(
                tstep,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cell_state.group_ix as usize],
                rng,
            )?;
            // println!(
            //     "                      old_v = {},\n                      new_v = {}, \n                      delta_mag: {}",
            //     cell_state.core.vertex_coords[vi],
            //     new_cell_state.core.vertex_coords[vi],
            //     (cell_state.core.vertex_coords[vi]
            //         - new_cell_state.core.vertex_coords[vi])
            //         .mag()
            // );

            // println!("----------------------------");

            // if (tstep == 347 && ci == 1) || (tstep == 346 && ci == 0)
            // {
            //     println!(
            //         "cell{} at {}->{}: {}",
            //         ci,
            //         tstep,
            //         tstep + 1,
            //         poly_to_string(
            //             &new_cell_state.state.vertex_coords
            //         )
            //     );
            // }
            // println!("----------------------");
            interaction_generator
                .update(ci, &new_cell_state.core.poly);
            confirm_volume_exclusion(
                &new_cell_state.core.poly,
                &interaction_generator.get_contact_data(ci),
                &format!("world cell {}", ci),
            )?;
            new_cell_states[ci] = new_cell_state;
        }
        Ok(Cells {
            cell_states: new_cell_states,
            interactions: interaction_generator.generate(),
        })
    }
}

#[derive(Deserialize, Serialize)]
pub struct WorldHistory {
    tstep: u32,
    interaction_generator: InteractionGenerator,
    rng: Pcg32,
    cells: Cells,
}

pub struct World {
    tstep_length: f32,
    tstep: u32,
    world_params: WorldParameters,
    group_params: Vec<Parameters>,
    pub history: Vec<WorldHistory>,
    cells: Cells,
    interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    output_dir: PathBuf,
    title: String,
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
    pub fn new(experiment: Experiment, output_dir: PathBuf) -> World {
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
            cell_group_ixs
                .append(&mut vec![gix; cg.num_cells as usize]);
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
        // Create initial Rho GTPase distributions.
        let cell_rgtps = cell_group_ixs
            .iter()
            .zip(cell_core_states.iter())
            .map(|(&gix, state)| {
                let parameters = &group_params[gix];
                state.calc_crl_rgtp_state(parameters)
            })
            .collect::<Vec<[RgtpActivityDiff; NVERTS]>>();
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
            cell_states.push(CellState::new(
                cell_ix as u32,
                group_ix as u32,
                cell_core_states[cell_ix],
                &cell_interactions[cell_ix],
                parameters,
                &mut cell_rng,
            ));
        }
        let cells = Cells {
            cell_states,
            interactions: cell_interactions.clone(),
        };
        let history = vec![WorldHistory {
            tstep: 0,
            interaction_generator: interaction_generator.clone(),
            rng,
            cells: cells.clone(),
        }];
        World {
            tstep: 0,
            tstep_length: char_quants.time(),
            world_params,
            group_params,
            history,
            cells,
            interaction_generator,
            rng,
            output_dir,
            title: experiment.title,
        }
    }

    pub fn as_history(&self) -> WorldHistory {
        WorldHistory {
            tstep: self.tstep,
            interaction_generator: self.interaction_generator.clone(),
            rng: self.rng,
            cells: self.cells.clone(),
        }
    }

    pub fn simulate(&mut self, final_tpoint: f32) {
        let num_tsteps =
            (final_tpoint / self.tstep_length).ceil() as u32;
        while self.tstep < num_tsteps {
            // println!(
            //     "tstep: {}/{}",
            //     self.tstep, num_tsteps
            // );
            #[cfg(feature = "debug_mode")]
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
                    self.save_history();
                    panic!("tstep: {}\n{}", self.tstep, e);
                });

            #[cfg(not(feature = "debug_mode"))]
            let new_cells = self.cells.simulate(
                &mut self.rng,
                self.cell_rngs.as_mut_slice(),
                &self.world_params,
                &self.group_params,
                &mut self.interaction_generator,
            );
            self.cells = new_cells;
            if self.tstep % 10 == 0 {
                self.history.push(self.as_history());
            }
            self.tstep += 1;
        }
    }

    pub fn save_history(&self) {
        dbg!("saving history...");
        #[cfg(feature = "debug_mode")]
        let name = format!("history_dbg_{}", &self.title);
        #[cfg(not(feature = "debug_mode"))]
        let name = format!("history_{}", &self.title);
        let mut path = self.output_dir.clone();
        path.push(format!("{}.json", &name));

        let mut f = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(path)
            .unwrap();

        for wh in self.history.iter() {
            serde_json::to_writer(&mut f, wh).unwrap();
        }
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
