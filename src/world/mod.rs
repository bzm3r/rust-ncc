// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod hardio;
#[cfg(feature = "validation")]
use crate::cell::confirm_volume_exclusion;
use crate::cell::core_state::CoreState;
use crate::cell::CellState;
use crate::experiments::{CellGroup, Experiment};
use crate::interactions::{
    CellInteractions, InteractionGenerator, RelativeRgtpActivity,
};
use crate::math::v2d::V2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::NVERTS;
//use rand_core::SeedableRng;
use crate::utils::pcg32::Pcg32;
use crate::world::hardio::{save_compact, save_full, Format};
use rand::seq::SliceRandom;
use rand::{RngCore, SeedableRng};
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::f64::consts::PI;
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize)]
pub struct Cells {
    pub cell_states: Vec<CellState>,
    pub interactions: Vec<CellInteractions>,
}

impl Cells {
    fn simulate(
        &self,
        tstep: usize,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
    ) -> Result<Cells, String> {
        let mut new_cell_states =
            vec![self.cell_states[0]; self.cell_states.len()];
        let shuffled_cells = {
            let mut crs =
                self.cell_states.iter().collect::<Vec<&CellState>>();
            crs.shuffle(rng);
            crs
        };
        for cell_state in shuffled_cells {
            let ci = cell_state.ix as usize;
            let contact_data =
                interaction_generator.get_contact_data(ci);

            let new_cell_state = cell_state.simulate_rkdp5(
                tstep,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cell_state.group_ix as usize],
                rng,
            )?;

            interaction_generator.update(
                tstep,
                ci,
                &new_cell_state.core.poly,
            );

            #[cfg(feature = "validation")]
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

#[derive(Deserialize, Serialize, Clone)]
pub struct Snapshot {
    pub tstep: usize,
    pub interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    pub cells: Cells,
}

pub struct World {
    tstep_length: f64,
    tstep: usize,
    world_params: WorldParameters,
    group_params: Vec<Parameters>,
    pub history: Vec<Snapshot>,
    cells: Cells,
    interaction_generator: InteractionGenerator,
    pub rng: Pcg32,
    out_dir: PathBuf,
    file_name: String,
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
            cell_states.push(CellState::new(
                cell_ix as usize,
                group_ix as usize,
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
        let history = vec![Snapshot {
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
            out_dir: output_dir,
            file_name: experiment.file_name,
        }
    }

    pub fn as_history(&self) -> Snapshot {
        Snapshot {
            tstep: self.tstep,
            interaction_generator: self.interaction_generator.clone(),
            rng: self.rng,
            cells: self.cells.clone(),
        }
    }

    pub fn simulate(
        &mut self,
        final_tpoint: f64,
        save_frequency: usize,
    ) {
        let num_tsteps =
            (final_tpoint / self.tstep_length).ceil() as usize;
        while self.tstep < num_tsteps {
            //println!("------------------");
            let tstep = self.tstep;
            //println!("tstep: {}", self.tstep);
            // if self.tstep == 1760 {
            //     println!("at tstep 1760!")
            // }
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
            if self.tstep % save_frequency == 0 {
                self.history.push(self.as_history());
            }
            // if true {
            //     let ics = self.interaction_generator.generate();
            //     let c0_cps = self
            //         .interaction_generator
            //         .get_close_points(0, 4, 1);
            //     let c1_cps = self
            //         .interaction_generator
            //         .get_close_points(1, 12, 0);
            //
            //     if (c0_cps.len() > 0 || c1_cps.len() > 0)
            //         && (c0_cps
            //             .iter()
            //             .any(|cp| cp.get_vector_to_mag() < 0.5)
            //             || c1_cps
            //                 .iter()
            //                 .any(|cp| cp.get_vector_to_mag() < 0.5))
            //     {
            //         println!("tstep = {}", tstep);
            //         println!(
            //             "ci (0, 4, oci: 1) = {}\nx_adhs[4] = {}",
            //             c0_cps
            //                 .iter()
            //                 .map(|cp| format!("{}", cp))
            //                 .collect::<Vec<String>>()
            //                 .join(", "),
            //             ics[0].x_adhs[4]
            //         );
            //         println!(
            //             "ci (1, 12, oci: 0) = {}\nx_adh[12] = {}",
            //             c1_cps
            //                 .iter()
            //                 .map(|cp| format!("{}", cp))
            //                 .collect::<Vec<String>>()
            //                 .join(", "),
            //             ics[1].x_adhs[12]
            //         );
            //         println!("------------------")
            //     }
            // }
        }
    }

    pub fn save_history(
        &self,
        compact: bool,
        formats: Vec<Format>,
    ) -> Result<(), Box<dyn Error>> {
        if compact {
            let data = self
                .history
                .iter()
                .map(|h| (h.tstep, h.cells.clone()))
                .collect::<Vec<(usize, Cells)>>();
            save_compact(
                data,
                &self.out_dir,
                formats,
                &self.file_name,
            )?;
        } else {
            let data = self
                .history
                .iter()
                .enumerate()
                .map(|(ix, snapshot)| (ix as usize, snapshot.clone()))
                .collect::<Vec<(usize, Snapshot)>>();
            save_full(data, &self.out_dir, formats, &self.file_name)?;
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
                + (row as f64) * row_delta
                + (col as f64) * col_delta;
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
