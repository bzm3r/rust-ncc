// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod hardio;
pub mod interactions;

use crate::cell::chemistry::RacRandState;
use crate::cell::state::{ChemState, GeomState, MechState, State};
use crate::cell::{Cell, VertexGenInfo};
use crate::experiment::{Experiment, GroupLayout};
use crate::math::p2d::P2D;
use crate::parameters::quantity::Length;
use crate::parameters::{InputParameters, Parameters, WorldParameters};
use crate::world::hardio::{save_data, save_schema};
use crate::world::interactions::InteractionState;
use crate::NVERTS;
use avro_rs::Writer;
use avro_schema_derive::Schematize;
use once_cell::unsync::Lazy;
use rand::rngs::SmallRng;
use rand::{thread_rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldState {
    tstep: u32,
    cells: Vec<Cell>,
    interactions: Vec<InteractionState>,
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldMechState {
    tstep: u32,
    state: Vec<MechState>,
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldRacRandState {
    tstep: u32,
    state: Vec<RacRandState>,
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldGeomState {
    tstep: u32,
    state: Vec<GeomState>,
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldChemState {
    tstep: u32,
    state: Vec<ChemState>,
}

impl WorldState {
    fn simulate(
        &self,
        cell_rngs: &mut [Option<RandomEventGenerator>],
        group_parameters: &[Parameters],
    ) -> WorldState {
        let mut cells = Vec::with_capacity(self.cells.len());
        for c in self.cells.iter() {
            cells.push(c.simulate_euler(
                self.tstep,
                &self.interactions[c.ix as usize],
                cell_rngs[c.ix as usize].as_mut(),
                &group_parameters[c.group_ix as usize],
            ));
        }
        let interactions = WorldState::calc_interactions(&cells);
        WorldState {
            tstep: self.tstep + 1,
            cells,
            interactions,
        }
    }

    fn calc_interactions(cells: &[Cell]) -> Vec<InteractionState> {
        //let contact_mat = ContactMatrix::calc(&cells);
        let mut r = vec![];
        for _ in cells.iter() {
            let x_cil = [0.0; NVERTS as usize];
            r.push(InteractionState {
                x_cil,
                x_chemoas: [0.0; NVERTS as usize],
                x_coas: [0.0; NVERTS as usize],
                x_bdrys: [0.0; NVERTS as usize],
            })
        }
        r
    }

    pub fn extract_geom_state(&self) -> WorldGeomState {
        WorldGeomState {
            tstep: self.tstep,
            state: self
                .cells
                .iter()
                .map(|c| State::calc_geom_state(&c.state))
                .collect::<Vec<GeomState>>(),
        }
    }

    pub fn extract_mech_state(&self, group_parameters: &[Parameters]) -> WorldMechState {
        WorldMechState {
            tstep: self.tstep,
            state: self
                .cells
                .iter()
                .map(|c| {
                    let gs = State::calc_geom_state(&c.state);
                    State::calc_mech_state(&c.state, &gs, &group_parameters[c.group_ix as usize])
                })
                .collect::<Vec<MechState>>(),
        }
    }

    pub fn extract_chem_state(&self, group_parameters: &[Parameters]) -> WorldChemState {
        WorldChemState {
            tstep: self.tstep,
            state: self
                .cells
                .iter()
                .map(|c| {
                    let gs = State::calc_geom_state(&c.state);
                    let ms = State::calc_mech_state(
                        &c.state,
                        &gs,
                        &group_parameters[c.group_ix as usize],
                    );
                    State::calc_chem_state(
                        &c.state,
                        &gs,
                        &ms,
                        &c.rac_rand_state,
                        &self.interactions[c.ix as usize],
                        &group_parameters[c.group_ix as usize],
                    )
                })
                .collect::<Vec<ChemState>>(),
        }
    }

    #[allow(unused)]
    pub fn extract_rac_rand_state(&self) -> WorldRacRandState {
        WorldRacRandState {
            tstep: self.tstep,
            state: self
                .cells
                .iter()
                .map(|c| c.rac_rand_state)
                .collect::<Vec<RacRandState>>(),
        }
    }
}

pub struct RandomEventGenerator {
    pub(crate) rng: SmallRng,
    distrib: Normal<f32>,
}

impl RandomEventGenerator {
    pub fn sample(&mut self) -> f32 {
        self.distrib.sample(&mut self.rng)
    }
}

pub struct World {
    world_parameters: WorldParameters,
    group_parameters: Vec<Parameters>,
    history: Vec<WorldState>,
    cell_regs: Vec<Option<RandomEventGenerator>>,
    state: WorldState,
}

impl World {
    pub fn new(experiment: Experiment) -> World {
        let mut num_cells = 0;
        let mut cells = vec![];
        let world_parameters = WorldParameters::default();
        let mut group_parameters = vec![];
        let mut primary = thread_rng();
        let mut cell_regs = vec![];
        for (gix, cg) in experiment.cell_groups.iter().enumerate() {
            let user_parameters =
                InputParameters::default().apply_overrides(&cg.parameter_overrides);
            let parameters = user_parameters.gen_parameters(&world_parameters);
            let cell_centroids = gen_cell_centroids(
                &cg.layout,
                cg.num_cells,
                user_parameters.cell_diam,
                &world_parameters,
            )
            .unwrap();
            for cc in cell_centroids.into_iter() {
                let normal = Lazy::new(|| {
                    Normal::new(parameters.rand_avg_t, parameters.rand_std_t).unwrap()
                });
                let mut creg = if parameters.randomization {
                    Some(RandomEventGenerator {
                        rng: SmallRng::from_rng(&mut primary).unwrap(),
                        distrib: *normal,
                    })
                } else {
                    None
                };
                cells.push(Cell::new(
                    num_cells,
                    gix as u32,
                    VertexGenInfo::Centroid(cc),
                    &parameters,
                    creg.as_mut(),
                ));
                cell_regs.push(creg);
                num_cells += 1;
            }
            group_parameters.push(parameters);
        }
        let interactions = WorldState::calc_interactions(&cells);
        let state = WorldState {
            tstep: 0,
            cells,
            interactions,
        };
        let history = vec![state.clone()];
        World {
            world_parameters,
            group_parameters,
            cell_regs,
            history,
            state,
        }
    }

    pub fn simulate(&mut self, final_tpoint: f32) {
        let num_tsteps = (final_tpoint / self.world_parameters.time()).ceil() as u32;

        while self.state.tstep < num_tsteps {
            // println!("========================================");
            // println!("tstep: {}/{}", self.state.tstep, num_tsteps);
            let new_state = self
                .state
                .simulate(self.cell_regs.as_mut_slice(), &self.group_parameters);
            self.history.push(new_state.clone());
            self.state = new_state;
            //println!("========================================")
        }
    }

    pub fn save_history(&self, output_dir: &PathBuf) {
        let name = "history";
        let schema = WorldState::schematize(None);
        save_schema(name, &schema, output_dir);
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer.append_ser(ws).unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        save_data(name, &encoded, output_dir);
    }

    pub fn save_mech_history(&self, output_dir: &PathBuf) {
        let name = "mech_hist";
        let schema = WorldMechState::schematize(None);
        save_schema(name, &schema, output_dir);
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer
                .append_ser(ws.extract_mech_state(&self.group_parameters))
                .unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        save_data(name, &encoded, output_dir);
    }

    #[allow(unused)]
    pub fn save_chem_history(&self, output_dir: &PathBuf) {
        let name = "chem_hist";
        let schema = WorldChemState::schematize(None);
        save_schema(name, &schema, output_dir);
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer
                .append_ser(ws.extract_chem_state(&self.group_parameters))
                .unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        save_data(name, &encoded, output_dir);
    }

    pub fn save_geom_history(&self, output_dir: &PathBuf) {
        let name = "geom_hist";
        let schema = WorldGeomState::schematize(None);
        save_schema(name, &schema, output_dir);
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer.append_ser(ws.extract_geom_state()).unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        save_data(name, &encoded, output_dir);
    }

    #[allow(unused)]
    pub fn save_rand_history(&self, output_dir: &PathBuf) {
        let name = "rand_hist";
        let schema = WorldRacRandState::schematize(None);
        save_schema(name, &schema, output_dir);
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer
                .append_ser(ws.extract_chem_state(&self.group_parameters))
                .unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        save_data(name, &encoded, output_dir);
    }
}

fn gen_cell_centroids(
    layout: &GroupLayout,
    num_cells: u32,
    cell_diam: Length,
    world_parameters: &WorldParameters,
) -> Result<Vec<P2D>, String> {
    let cell_diam = world_parameters.normalize(&cell_diam);
    let group_centroid = P2D {
        x: world_parameters.normalize(&layout.centroid[0]),
        y: world_parameters.normalize(&layout.centroid[1]),
    };
    if layout.width * layout.height >= num_cells {
        let nr = layout.height / num_cells;
        let mut r = vec![];
        let first_cell_centroid = {
            let delta = P2D {
                x: ((layout.width - 1) as f32) * 0.5 * cell_diam,
                y: ((layout.height - 1) as f32) * 0.5 * cell_diam,
            };
            group_centroid - delta
        };
        let cd = P2D {
            x: 0.5 * cell_diam,
            y: 0.5 * cell_diam,
        };
        for ix in 0..num_cells {
            let row = (ix / nr) as f32;
            let col = (ix as f32) - row;
            r.push(first_cell_centroid + cd.scalar_mul_x(row) + cd.scalar_mul_y(col));
        }
        Ok(r)
    } else {
        Err(format!(
            "Cell group layout area ({}x{}) not large enough to fit {} cells.",
            layout.width, layout.height, num_cells
        ))
    }
}
