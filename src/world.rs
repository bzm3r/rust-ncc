// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::cell::Cell;
use crate::experiment::{Experiment, GroupLayout};
use crate::math::P2D;
use crate::parameters::{InputParameters, Parameters, WorldParameters};
use crate::quantity::Length;
use avro_rs::{Schema, Writer};
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::fs::{OpenOptions};
use std::io::Write;
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldState {
    tstep: u32,
    cells: Vec<Cell>,
}

impl WorldState {
    fn simulate(&self, group_parameters: &[Parameters]) -> WorldState {
        let mut new_cells = Vec::with_capacity(self.cells.len());
        for cell in self.cells.iter() {
            let new_cell = cell.simulate_rkdp5(self.tstep, &group_parameters[cell.group_ix as usize]);
            new_cells.push(new_cell);
        }
        WorldState {
            tstep: self.tstep + 1,
            cells: new_cells,
        }
    }
}

pub struct World {
    world_parameters: WorldParameters,
    group_parameters: Vec<Parameters>,
    history: Vec<WorldState>,
    state: WorldState,
    state_schema: Schema,
    output_dir: PathBuf,
}

impl World {
    pub fn new(experiment: Experiment, output_dir: PathBuf) -> World {
        let mut num_cells = 0;
        let mut cells = vec![];
        let world_parameters = WorldParameters::default();
        let mut group_parameters = vec![];
        for (gix, cg) in experiment.cell_groups.iter().enumerate() {
            let user_params = InputParameters::default().apply_overrides(&cg.parameter_overrides);
            let params = user_params.gen_params(&world_parameters);
            let cell_centroids = gen_cell_centroids(
                &cg.layout,
                cg.num_cells,
                user_params.cell_diam,
                &world_parameters,
            )
            .unwrap();
            for cc in cell_centroids.into_iter() {
                cells.push(Cell::new(num_cells, gix as u32, &params, cc));
                num_cells += 1;
            }
            group_parameters.push(params);
        }

        let state = WorldState { tstep: 0, cells };

        let history = vec![state.clone()];
        let raw_schema = WorldState::raw_schema();
        {
            let mut avsc_path = output_dir.clone();
            avsc_path.push("schema");
            avsc_path.set_extension("avsc");

            let mut f = OpenOptions::new().write(true).create(true).truncate(true).open(avsc_path).unwrap();
            f.write_all(raw_schema.as_bytes()).unwrap();
        }
        let state_schema = Schema::parse_str(&raw_schema).unwrap();

        World {
            world_parameters,
            group_parameters,
            history,
            state,
            state_schema,
            output_dir,
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

    pub fn save_history(&self) {
        let mut writer = Writer::new(&self.state_schema, Vec::new());//Writer::with_codec(&self.state_schema, Vec::new(), avro_rs::Codec::Deflate);
        for ws in self.history.iter() {
            writer.append_ser(ws).unwrap();
        }
        let encoded = writer.into_inner().unwrap();
        {
            let mut avro_path = self.output_dir.clone();
            avro_path.push("dat");
            avro_path.set_extension("avro");

            let mut f = OpenOptions::new().write(true).create(true).truncate(true).open(avro_path).unwrap();
            f.write_all(&encoded).unwrap();
        }
    }
}

fn gen_cell_centroids(
    layout: &GroupLayout,
    num_cells: u32,
    cell_diam: Length,
    world_params: &WorldParameters,
) -> Result<Vec<P2D>, String> {
    let cell_diam = world_params.normalize(&cell_diam);
    let group_centroid = P2D {
        x: world_params.normalize(&layout.centroid[0]),
        y: world_params.normalize(&layout.centroid[1]),
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
