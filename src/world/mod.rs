// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod hardio;
use crate::cell::ModelCell;
use crate::experiments::{CellGroup, Experiment};
use crate::interactions::{find_contacts, Contacts, InteractionState};
use crate::math::geometry::BBox;
use crate::math::p2d::P2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::world::hardio::{save_data, save_schema};
use crate::NVERTS;
use avro_rs::Writer;
use avro_schema_derive::Schematize;
use once_cell::unsync::Lazy;
use rand::rngs::SmallRng;
use rand::{thread_rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;
use std::path::PathBuf;

#[derive(Clone, Deserialize, Serialize, Schematize)]
struct WorldState {
    tstep: u32,
    cells: Vec<ModelCell>,
    interactions: Vec<InteractionState>,
    #[serde(skip)]
    contacts: Contacts,
}

impl WorldState {
    fn simulate(
        &self,
        cell_regs: &mut [Option<RandomEventGenerator>],
        world_parameters: &WorldParameters,
        group_parameters: &[Parameters],
    ) -> WorldState {
        let cell_vcs = self
            .cells
            .iter()
            .map(|c| c.state.vertex_coords)
            .collect::<Vec<[P2D; NVERTS]>>();
        let cell_bbs = cell_vcs
            .iter()
            .map(|vcs| BBox::calc(vcs))
            .collect::<Vec<BBox>>();
        let cell_polys = cell_bbs
            .into_iter()
            .zip(cell_vcs.into_iter())
            .collect::<Vec<(BBox, [P2D; NVERTS])>>();
        let mut cells = Vec::with_capacity(self.cells.len());
        for (ci, c) in self.cells.iter().enumerate() {
            //println!("cell: {}", ix);
            let contact_polys = self
                .contacts
                .get_contacts(ci)
                .into_iter()
                .map(|ci| &cell_polys[ci])
                .collect::<Vec<&(BBox, [P2D; NVERTS])>>();
            cells.push(c.simulate_euler(
                self.tstep,
                &self.interactions[ci],
                contact_polys,
                cell_regs[ci].as_mut(),
                world_parameters,
                &group_parameters[c.group_ix as usize],
            ));
        }
        let cell_vcs = cells
            .iter()
            .map(|c| c.state.vertex_coords)
            .collect::<Vec<[P2D; NVERTS]>>();
        let cell_bboxes = cell_vcs
            .iter()
            .map(|vcs| BBox::calc(vcs))
            .collect::<Vec<BBox>>();
        let contacts = find_contacts(&cell_bboxes, world_parameters.close_criterion);
        let contact_dists = contacts.calc_dists(cell_vcs.as_slice());
        WorldState {
            tstep: self.tstep + 1,
            cells,
            interactions: contact_dists.interactions(&world_parameters.cil),
            contacts,
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
    tstep_length: f32,
    world_parameters: WorldParameters,
    group_parameters: Vec<Parameters>,
    history: Vec<WorldState>,
    cell_regs: Vec<Option<RandomEventGenerator>>,
    state: WorldState,
}

fn gen_vertex_coords(centroid: &P2D, radius: f32) -> [P2D; NVERTS] {
    let mut r = [P2D::default(); NVERTS];
    (0..NVERTS).for_each(|vix| {
        let vf = (vix as f32) / (NVERTS as f32);
        let theta = 2.0 * PI * vf;
        r[vix] = P2D {
            x: centroid.x + theta.cos() * radius,
            y: centroid.y + theta.sin() * radius,
        };
    });
    r
}

impl World {
    pub fn new(experiment: Experiment) -> World {
        let Experiment {
            basic_quants,
            world_parameters,
            cell_groups,
            ..
        } = experiment;
        let group_parameters = cell_groups
            .iter()
            .map(|cg| cg.parameters.clone())
            .collect::<Vec<Parameters>>();

        let mut cell_group_ixs = vec![];
        let mut cell_centroids = vec![];
        cell_groups.iter().enumerate().for_each(|(gix, cg)| {
            cell_group_ixs.append(&mut vec![gix; cg.num_cells as usize]);
            cell_centroids.append(&mut gen_cell_centroids(cg).unwrap())
        });
        let cell_vcs = cell_group_ixs
            .iter()
            .zip(cell_centroids.iter())
            .map(|(&gix, cc)| gen_vertex_coords(cc, group_parameters[gix].cell_r))
            .collect::<Vec<[P2D; NVERTS]>>();
        let cell_bboxes = cell_vcs
            .iter()
            .map(|vcs| BBox::calc(vcs))
            .collect::<Vec<BBox>>();
        let contacts = find_contacts(&cell_bboxes, world_parameters.close_criterion);
        let contact_dists = contacts.calc_dists(&cell_vcs);
        let interactions = contact_dists.interactions(&world_parameters.cil);
        let mut cells = vec![];
        let mut primary = thread_rng();
        let mut cell_regs = vec![];
        for (ix, gix) in cell_group_ixs.into_iter().enumerate() {
            let parameters = &group_parameters[gix];
            let normal =
                Lazy::new(|| Normal::new(parameters.rand_avg_t, parameters.rand_std_t).unwrap());
            let mut creg = if parameters.randomization {
                Some(RandomEventGenerator {
                    rng: SmallRng::from_rng(&mut primary).unwrap(),
                    distrib: *normal,
                })
            } else {
                None
            };
            cells.push(ModelCell::new(
                ix as u32,
                gix as u32,
                cell_vcs[ix],
                &interactions[ix],
                parameters,
                creg.as_mut(),
            ));
            cell_regs.push(creg);
        }
        let state = WorldState {
            tstep: 0,
            cells,
            interactions,
            contacts,
        };
        let history = vec![state.clone()];
        World {
            tstep_length: basic_quants.time(),
            world_parameters,
            group_parameters,
            cell_regs,
            history,
            state,
        }
    }

    pub fn simulate(&mut self, final_tpoint: f32) {
        let num_tsteps = (final_tpoint / self.tstep_length).ceil() as u32;
        while self.state.tstep < num_tsteps {
            // println!("========================================");
            // println!("tstep: {}/{}", self.state.tstep, num_tsteps);
            let new_state = self.state.simulate(
                self.cell_regs.as_mut_slice(),
                &self.world_parameters,
                &self.group_parameters,
            );
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
}

fn gen_cell_centroids(cg: &CellGroup) -> Result<Vec<P2D>, String> {
    let CellGroup {
        num_cells,
        layout,
        parameters,
    } = cg;
    let cell_r = parameters.cell_r;
    if layout.width * layout.height >= *num_cells {
        let mut r = vec![];
        let first_cell_centroid = P2D {
            x: layout.bottom_left.x + cell_r,
            y: layout.bottom_left.y + cell_r,
        };
        let delta_x = P2D {
            x: 2.0 * cell_r,
            y: 0.0,
        };
        let delta_y = P2D {
            x: 0.0,
            y: 2.0 * cell_r,
        };
        for ix in 0..*num_cells {
            let row = (ix / layout.width) as f32;
            let col = ix as f32 - row;
            r.push(first_cell_centroid + col * delta_x + row * delta_y);
        }
        Ok(r)
    } else {
        Err(format!(
            "Cell group layout area ({}x{}) not large enough to fit {} cells.",
            layout.width, layout.height, num_cells
        ))
    }
}
