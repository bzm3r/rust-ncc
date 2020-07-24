// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::cell::{Cell, CellState, ChemState, GeomState, MechState, VertexGenInfo};
use crate::consts::NVERTS;
use crate::experiment::{Experiment, GroupLayout};
use crate::math::P2D;
use crate::parameters::{InputParameters, Parameters, WorldParameters};
use crate::quantity::Length;
use avro_rs::{Schema, Writer};
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::io::Write;
use std::path::PathBuf;

#[derive(Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize)]
pub struct InteractionState {
    pub(crate) x_cils: [f32; NVERTS as usize],
    pub(crate) x_chemoas: [f32; NVERTS as usize],
    pub(crate) x_coas: [f32; NVERTS as usize],
    x_bdrys: [f32; NVERTS as usize],
}

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
struct WorldGeomState {
    tstep: u32,
    state: Vec<GeomState>,
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldChemState {
    tstep: u32,
    chem_state: Vec<ChemState>,
}

impl WorldState {
    fn simulate(&self, group_parameters: &[Parameters]) -> WorldState {
        let mut cells = Vec::with_capacity(self.cells.len());
        for c in self.cells.iter() {
            cells.push(c.simulate_rkdp5(
                self.tstep,
                &self.interactions[c.ix as usize],
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
        let mut r = vec![];
        for _ in cells.iter() {
            r.push(InteractionState::default())
        }
        r
    }

    pub fn extract_geom_state(&self) -> WorldGeomState {
        WorldGeomState {
            tstep: self.tstep,
            state: self
                .cells
                .iter()
                .map(|c| CellState::calc_geom_state(&c.state))
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
                    let gs = CellState::calc_geom_state(&c.state);
                    CellState::calc_mech_state(
                        &c.state,
                        &gs,
                        &group_parameters[c.group_ix as usize],
                    )
                })
                .collect::<Vec<MechState>>(),
        }
    }

    pub fn extract_chem_state(&self, group_parameters: &[Parameters]) -> WorldChemState {
        WorldChemState {
            tstep: self.tstep,
            chem_state: self
                .cells
                .iter()
                .map(|c| {
                    let gs = CellState::calc_geom_state(&c.state);
                    CellState::calc_chem_state(
                        &c.state,
                        &gs,
                        &c.rac_randomization,
                        &self.interactions[c.ix as usize],
                        &group_parameters[c.group_ix as usize],
                    )
                })
                .collect::<Vec<ChemState>>(),
        }
    }
}

pub struct World {
    world_parameters: WorldParameters,
    group_parameters: Vec<Parameters>,
    history: Vec<WorldState>,
    state: WorldState,
}

impl World {
    pub fn new(experiment: Experiment) -> World {
        let mut num_cells = 0;
        let mut cells = vec![];
        let world_parameters = WorldParameters::default();
        let mut group_parameters = vec![];
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
                cells.push(Cell::new(
                    num_cells,
                    gix as u32,
                    VertexGenInfo::Centroid(cc),
                    &parameters,
                ));
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
            history,
            state,
        }
    }

    pub fn simulate(&mut self, final_tpoint: f32) {
        let num_tsteps = (final_tpoint / self.world_parameters.time()).ceil() as u32;

        while self.state.tstep < num_tsteps {
            let new_state = self.state.simulate(&self.group_parameters);
            self.history.push(new_state.clone());
            self.state = new_state;
        }
    }

    pub fn save_history(&self, output_dir: &PathBuf) {
        let name = "history";
        let raw_schema = WorldState::raw_schema();
        World::save_schema(name, &raw_schema, output_dir);
        let schema = Schema::parse_str(&raw_schema).unwrap();
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer.append_ser(ws).unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        World::save_data(name, &encoded, output_dir);
    }

    pub fn save_mech_history(&self, output_dir: &PathBuf) {
        let name = "mech_hist";
        let raw_schema = WorldMechState::raw_schema();
        World::save_schema(name, &raw_schema, output_dir);
        let schema = Schema::parse_str(&raw_schema).unwrap();
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer
                .append_ser(ws.extract_mech_state(&self.group_parameters))
                .unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        World::save_data(name, &encoded, output_dir);
    }

    pub fn save_chem_history(&self, output_dir: &PathBuf) {
        let name = "chem_hist";
        let raw_schema = WorldChemState::raw_schema();
        World::save_schema(name, &raw_schema, output_dir);
        let schema = Schema::parse_str(&raw_schema).unwrap();
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer
                .append_ser(ws.extract_chem_state(&self.group_parameters))
                .unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        World::save_data(name, &encoded, output_dir);
    }

    pub fn save_geom_history(&self, output_dir: &PathBuf) {
        let name = "geom_hist";
        let raw_schema = WorldGeomState::raw_schema();
        World::save_schema(name, &raw_schema, output_dir);
        let schema = Schema::parse_str(&raw_schema).unwrap();
        let mut writer = Writer::new(&schema, Vec::new()); //Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer.append_ser(ws.extract_geom_state()).unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        World::save_data(name, &encoded, output_dir);
    }

    pub fn save_schema(name: &str, raw_schema: &str, output_dir: &PathBuf) {
        let mut avsc_path = output_dir.clone();
        avsc_path.push(format!("{}_schema", name));
        avsc_path.set_extension("avsc");

        let mut f = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(avsc_path)
            .unwrap();
        f.write_all(raw_schema.as_bytes()).unwrap();
    }

    pub fn save_data(name: &str, encoded: &[u8], output_dir: &PathBuf) {
        let mut path = output_dir.clone();
        path.push(format!("{}_dat", name));
        path.set_extension("avro");

        let mut f = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(path)
            .unwrap();
        f.write_all(&encoded).unwrap();
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
