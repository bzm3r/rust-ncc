pub mod py_comp;

// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::cell::states::Core;
use crate::cell::Cell;
use crate::exp_setup::{CellGroup, Experiment};
use crate::hardio::AsyncWriter;
use crate::interactions::{
    InteractionGenerator, Interactions, RelativeRgtpActivity,
};
use crate::math::v2d::V2d;
use crate::parameters::quantity::Quantity;
use crate::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::world::py_comp::execute_py_model;
use crate::NVERTS;
use rand::seq::SliceRandom;
use rand::{RngCore, SeedableRng};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::path::PathBuf;

#[derive(
    Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct WorldCells {
    pub tpoint: f64,
    pub cells: Vec<Cell>,
    pub interactions: Vec<Interactions>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub struct EulerOpts {
    pub num_int_steps: usize,
}

impl EulerOpts {
    pub fn dt(&self) -> f64 {
        1.0 / (self.num_int_steps as f64)
    }
}

impl Default for EulerOpts {
    fn default() -> Self {
        EulerOpts { num_int_steps: 10 }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub struct Rkdp5Opts {
    pub max_iters: usize,
    pub atol: f64,
    pub rtol: f64,
    pub init_h_scale: f64,
}

impl Default for Rkdp5Opts {
    fn default() -> Self {
        Rkdp5Opts {
            max_iters: 20,
            atol: 1e-3,
            rtol: 1e-3,
            init_h_scale: 0.1,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
pub enum IntegratorOpts {
    Euler(EulerOpts),
    EulerDebug(EulerOpts),
    Rkdp5(Rkdp5Opts),
}

impl Default for IntegratorOpts {
    fn default() -> Self {
        IntegratorOpts::Rkdp5(Rkdp5Opts::default())
    }
}

impl WorldCells {
    fn simulate_rkdp5(
        &self,
        tpoint: f64,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
        int_opts: Rkdp5Opts,
    ) -> Result<WorldCells, String> {
        let mut new_cells = self.cells.clone();
        let shuffled_cells = {
            let mut crs = self.cells.iter().collect::<Vec<&Cell>>();
            crs.shuffle(rng);
            crs
        };
        let dt = 1.0;
        for cells in shuffled_cells {
            let ci = cells.ix;
            let contact_data =
                interaction_generator.get_contact_data(ci);

            let new_cell = cells.simulate_rkdp5(
                tpoint,
                dt,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cells.group_ix],
                rng,
                int_opts,
            )?;

            interaction_generator.update(ci, &new_cell.core.poly);

            new_cells[ci] = new_cell;
        }
        let rel_rgtps = new_cells
            .iter()
            .map(|c| {
                c.core.calc_relative_rgtp_activity(
                    &group_parameters[c.group_ix],
                )
            })
            .collect::<Vec<[RelativeRgtpActivity; NVERTS]>>();
        Ok(WorldCells {
            tpoint: tpoint + dt,
            cells: new_cells,
            interactions: interaction_generator.generate(&rel_rgtps),
        })
    }

    fn simulate_euler(
        &self,
        tpoint: f64,
        rng: &mut Pcg32,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
        int_opts: EulerOpts,
    ) -> Result<WorldCells, String> {
        let mut new_cells = self.cells.clone();
        let shuffled_cells = {
            let mut crs = self.cells.iter().collect::<Vec<&Cell>>();
            crs.shuffle(rng);
            crs
        };
        for cell in shuffled_cells {
            let ci = cell.ix;
            let contact_data =
                interaction_generator.get_contact_data(ci);

            let new_cell = cell.simulate_euler(
                tpoint,
                &self.interactions[ci],
                contact_data,
                world_parameters,
                &group_parameters[cell.group_ix],
                rng,
                int_opts,
            )?;

            interaction_generator.update(ci, &new_cell.core.poly);

            new_cells[ci] = new_cell;
        }
        let rel_rgtps = new_cells
            .iter()
            .map(|c| {
                c.core.calc_relative_rgtp_activity(
                    &group_parameters[c.group_ix],
                )
            })
            .collect::<Vec<[RelativeRgtpActivity; NVERTS]>>();
        Ok(WorldCells {
            tpoint: tpoint + 1.0,
            cells: new_cells,
            interactions: interaction_generator.generate(&rel_rgtps),
        })
    }

    fn simulate_euler_debug(
        &self,
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
        interaction_generator: &mut InteractionGenerator,
        int_opts: EulerOpts,
    ) -> Result<Vec<WorldCells>, String> {
        let mut out = vec![];
        let dt = int_opts.dt();
        for int_step in 0..(int_opts.num_int_steps + 1) {
            let mut dummy = self.clone();
            dummy.tpoint = self.tpoint + dt * int_step as f64;
            out.push(dummy);
        }
        let mut final_states = self.cells.clone();
        let last_ix = int_opts.num_int_steps - 1;
        for cell in self.cells.iter() {
            let rel_rgtps = final_states
                .iter()
                .map(|c| {
                    c.core.calc_relative_rgtp_activity(
                        &group_parameters[c.group_ix],
                    )
                })
                .collect::<Vec<[RelativeRgtpActivity; NVERTS]>>();
            let interactions =
                interaction_generator.generate(&rel_rgtps);
            let ci = cell.ix;
            let contact_data =
                interaction_generator.get_contact_data(ci);
            let this_interactions = &interactions[ci];

            let r = cell.simulate_euler_debug(
                this_interactions,
                contact_data,
                world_parameters,
                &group_parameters[cell.group_ix],
                int_opts,
            )?;
            for (int_step, int_state) in r.iter().enumerate() {
                out[int_step + 1].cells[ci] = *int_state;
                out[int_step + 1].interactions[ci] =
                    *this_interactions;
            }
            final_states[ci] = r[last_ix];
            interaction_generator
                .update(ci, &final_states[ci].core.poly);
        }
        Ok(out)
    }
}

#[derive(
    Deserialize, Serialize, Clone, Default, Debug, PartialEq,
)]
pub struct WorldInfo {
    pub final_t: f64,
    pub snap_period: f64,
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    pub cell_params: Vec<Parameters>,
}

impl Iterator for WorldCells {
    type Item = ();

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

#[derive(Clone)]
pub struct WorldState {
    pub tpoint: f64,
    pub cells: WorldCells,
    pub rng: Pcg32,
}

pub struct World {
    final_t: f64,
    char_quants: CharQuantities,
    state: WorldState,
    params: WorldParameters,
    cell_group_params: Vec<Parameters>,
    writer: Option<AsyncWriter>,
    interaction_generator: InteractionGenerator,
    snap_period: f64,
    int_opts: IntegratorOpts,
}

fn gen_poly(centroid: &V2d, radius: f64) -> [V2d; NVERTS] {
    let mut r = [V2d::default(); NVERTS];
    (0..NVERTS).for_each(|vix| {
        let vf = (vix as f64) / (NVERTS as f64);
        let theta = 2.0 * PI * vf;
        r[vix] = V2d {
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
            world_params,
            cell_groups,
            mut rng,
            final_t,
            snap_period,
            max_on_ram,
            int_opts,
            out_dir,
            py_main,
            name,
            ..
        } = experiment;
        let normed_final_t = char_quants.normalize(&final_t);
        let normed_snap_period = char_quants.normalize(&snap_period);
        let num_tsteps = normed_final_t.ceil() as usize;
        let expected_final_t = char_quants
            .normalize(&char_quants.t.scale(num_tsteps as f64));
        // Extract the parameters from each `CellGroup` object obtained
        // from the `Experiment`.
        let group_params = cell_groups
            .iter()
            .map(|cg| cg.parameters)
            .collect::<Vec<Parameters>>();

        // Create a list of indices of the groups. and create a vector
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
            .collect::<Vec<[V2d; NVERTS]>>();
        if let Some(pm) = &py_main {
            execute_py_model(
                &out_dir,
                &pm,
                &name,
                final_t.number(),
                snap_period.number(),
                cell_polys.len() as u32,
                world_params.interactions.phys_contact.cil_mag,
                world_params.interactions.coa.map(|coa_params| {
                    coa_params.vertex_mag * NVERTS as f64
                }),
            );
        }
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
        let mut cells = vec![];
        for (cell_ix, group_ix) in
            cell_group_ixs.into_iter().enumerate()
        {
            // Parameters that will be used by this cell. Determined
            // by figuring out which group it belongs to, as all cells
            // within a group use the same parameters.
            let parameters = &group_params[group_ix];
            let mut cell_rng = Pcg32::seed_from_u64(rng.next_u64());
            // Create a new cell.
            cells.push(Cell::new(
                cell_ix,
                group_ix,
                cell_core_states[cell_ix],
                parameters,
                &mut cell_rng,
            ));
        }
        let cells = WorldCells {
            tpoint: 0.0,
            cells,
            interactions: cell_interactions,
        };
        let writer = Some(Self::init_writer(
            out_dir,
            name,
            WorldInfo {
                final_t: expected_final_t,
                snap_period: normed_snap_period,
                char_quants,
                world_params: world_params.clone(),
                cell_params: cells
                    .cells
                    .iter()
                    .map(|s| group_params[s.group_ix])
                    .collect::<Vec<Parameters>>(),
            },
            max_on_ram,
        ));
        World {
            state: WorldState {
                tpoint: 0.0,
                cells,
                rng,
            },
            final_t: expected_final_t,
            char_quants,
            params: world_params,
            cell_group_params: group_params,
            interaction_generator,
            writer,
            int_opts,
            snap_period: normed_snap_period,
        }
    }

    pub fn save_state(&mut self) {
        if let Some(writer) = &mut self.writer {
            writer.push(self.state.clone());
        }
    }

    pub fn periodic_save(&mut self, last_saved: f64) -> f64 {
        if (self.state.tpoint - last_saved) > self.snap_period {
            if let Some(writer) = &mut self.writer {
                writer.push(self.state.clone());
            }
            self.state.tpoint
        } else {
            last_saved
        }
    }

    pub fn periodic_save_euler_debug(
        &mut self,
        curr_tpoint: f64,
        last_saved: f64,
    ) -> f64 {
        if (curr_tpoint - last_saved) >= self.snap_period {
            // println!(
            //     "curr_tpoint: {}, saving: {}",
            //     curr_tpoint, self.state.tpoint
            // );
            if let Some(writer) = &mut self.writer {
                writer.push(self.state.clone());
            }
            curr_tpoint
        } else {
            last_saved
        }
    }

    pub fn simulate_rkdp5(
        &mut self,
        save_cbor: bool,
        int_opts: Rkdp5Opts,
    ) {
        // Save initial state.
        self.save_state();
        let mut last_saved = self.state.tpoint;
        while self.state.tpoint < self.final_t {
            let new_cells: WorldCells = self
                .state
                .cells
                .simulate_rkdp5(
                    self.state.tpoint,
                    &mut self.state.rng,
                    &self.params,
                    &self.cell_group_params,
                    &mut self.interaction_generator,
                    int_opts,
                )
                .unwrap_or_else(|e| {
                    self.final_save(save_cbor, "panicking");
                    panic!("tstep: {}\n{}", self.state.tpoint, e);
                });

            self.state.tpoint = new_cells.tpoint;
            self.state.cells = new_cells;
            last_saved = self.periodic_save(last_saved);
        }
        self.final_save(save_cbor, "done");
    }

    pub fn simulate_euler(
        &mut self,
        save_cbor: bool,
        int_opts: EulerOpts,
    ) {
        // Save initial state.
        self.save_state();
        let mut last_saved = 0.0;
        while self.state.tpoint < self.final_t {
            let new_cells: WorldCells = self
                .state
                .cells
                .simulate_euler(
                    self.state.tpoint,
                    &mut self.state.rng,
                    &self.params,
                    &self.cell_group_params,
                    &mut self.interaction_generator,
                    int_opts,
                )
                .unwrap_or_else(|e| {
                    self.final_save(save_cbor, "panicking");
                    panic!("tstep: {}\n{}", self.state.tpoint, e);
                });

            self.state.tpoint = new_cells.tpoint;
            self.state.cells = new_cells;
            last_saved = self.periodic_save(last_saved);
        }
        self.final_save(save_cbor, "done");
    }

    pub fn simulate_euler_debug(
        &mut self,
        save_cbor: bool,
        int_opts: EulerOpts,
    ) {
        let mut last_saved = 0.0 - self.snap_period;
        while self.state.tpoint < self.final_t {
            let new_cells = self
                .state
                .cells
                .simulate_euler_debug(
                    &self.params,
                    &self.cell_group_params,
                    &mut self.interaction_generator,
                    int_opts,
                )
                .unwrap_or_else(|e| {
                    self.final_save(save_cbor, "panicking");
                    panic!("tstep: {}\n{}", self.state.tpoint, e);
                });
            let curr_tpoint = self.state.tpoint;
            let mut next_last_saved = 0.0;
            for cells in new_cells[..int_opts.num_int_steps].iter() {
                self.state.tpoint = cells.tpoint;
                self.state.cells = cells.clone();
                next_last_saved = self.periodic_save_euler_debug(
                    curr_tpoint,
                    last_saved,
                );
            }
            self.state.tpoint =
                new_cells[int_opts.num_int_steps].tpoint;
            self.state.cells =
                new_cells[int_opts.num_int_steps].clone();
            last_saved = next_last_saved;
        }
        self.final_save(save_cbor, "done");
    }

    pub fn simulate(&mut self, save_cbor: bool) {
        match self.int_opts {
            IntegratorOpts::Euler(int_opts) => {
                self.simulate_euler(save_cbor, int_opts)
            }
            IntegratorOpts::EulerDebug(int_opts) => {
                self.simulate_euler_debug(save_cbor, int_opts)
            }
            IntegratorOpts::Rkdp5(int_opts) => {
                self.simulate_rkdp5(save_cbor, int_opts)
            }
        }
    }

    pub fn info(&self) -> WorldInfo {
        WorldInfo {
            final_t: self.final_t,
            snap_period: self.snap_period,
            char_quants: self.char_quants,
            world_params: self.params.clone(),
            cell_params: self
                .state
                .cells
                .cells
                .iter()
                .map(|s| self.cell_group_params[s.group_ix])
                .collect::<Vec<Parameters>>(),
        }
    }

    pub fn init_writer(
        output_dir: PathBuf,
        file_name: String,
        info: WorldInfo,
        max_capacity: usize,
    ) -> AsyncWriter {
        AsyncWriter::new(
            output_dir,
            file_name,
            max_capacity,
            true,
            info,
        )
    }

    pub fn final_save(&mut self, save_cbor: bool, reason: &str) {
        if let Some(writer) = self.writer.take() {
            writer.finish(save_cbor, reason);
        }
    }
}

pub fn gen_cell_centroids(
    cg: &CellGroup,
) -> Result<Vec<V2d>, String> {
    let CellGroup {
        num_cells,
        layout,
        parameters,
    } = cg;
    let cell_r = parameters.cell_r;
    if layout.width * layout.height >= *num_cells {
        let mut r = vec![];
        let first_cell_centroid = V2d {
            x: layout.bottom_left.x + cell_r,
            y: layout.bottom_left.y + cell_r,
        };
        let row_delta = V2d {
            x: 0.0,
            y: 2.0 * cell_r,
        };
        let col_delta = V2d {
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
