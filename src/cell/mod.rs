// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod core_state;
pub mod mechanics;
pub mod rkdp5;

use crate::cell::chemistry::RacRandState;
use crate::cell::core_state::{ChemState, CoreState, DepStates, GeomState, MechState};
use crate::cell::rkdp5::AuxArgs;
use crate::interactions::{CloseCellInfo, InteractionState};
use crate::math::geometry::{calc_poly_area, move_inner_point_to_bdry, Bbox};
use crate::math::p2d::P2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::world::RandomEventGenerator;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

/// Model cell.
#[derive(Copy, Clone, Deserialize, Serialize, Schematize)]
pub struct ModelCell {
    /// Cell index.
    pub ix: u32,
    /// Index of cell type.
    pub group_ix: u32,
    pub state: CoreState,
    pub rac_rand_state: RacRandState,
    pub mech_state: MechState,
    pub geom_state: GeomState,
    pub chem_state: ChemState,
}

fn enforce_volume_exclusion(
    vcs: &[P2D; NVERTS],
    geom_state: &GeomState,
    close_cells: &CloseCellInfo,
    cell_bboxes: &[Bbox],
    cell_vcs: &[&[P2D; NVERTS]],
) -> [P2D; NVERTS] {
    let mut r = *vcs;
    for cix in 0..close_cells.num {
        let ovcs = cell_vcs[cix];
        let obbox = &cell_bboxes[cix];
        for (vix, vc) in vcs.iter().enumerate() {
            r[vix] = move_inner_point_to_bdry(vc, &geom_state.unit_inward_vecs[vix], obbox, ovcs);
        }
    }
    r
}

impl ModelCell {
    pub fn new(
        ix: u32,
        group_ix: u32,
        vertex_coords: [P2D; NVERTS],
        interactions: &InteractionState,
        parameters: &Parameters,
        reg: Option<&mut RandomEventGenerator>,
    ) -> ModelCell {
        let state = CoreState::new(vertex_coords, parameters.init_rac, parameters.init_rho);
        let rac_rand_state = RacRandState::init(
            match reg {
                Some(cr) => Some(&mut cr.rng),
                None => None,
            },
            parameters,
        );
        let DepStates {
            geom_state,
            chem_state,
            mech_state,
        } = CoreState::calc_dep_states(&state, &rac_rand_state, interactions, parameters);
        ModelCell {
            ix,
            group_ix,
            state,
            rac_rand_state,
            geom_state,
            chem_state,
            mech_state,
        }
    }

    #[allow(unused)]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        interactions: &InteractionState,
        close_cells: &CloseCellInfo,
        cell_bboxes: &[Bbox],
        cell_vcs: &[&[P2D; NVERTS]],
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> ModelCell {
        let mut state = self.state;
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for i in 0..nsteps {
            //println!("++++++++++++");
            //println!("{}", state);
            let dep_vars =
                CoreState::calc_dep_states(&state, &self.rac_rand_state, interactions, parameters);
            //println!("{}", dep_vars);
            //println!("++++++++++++");
            let delta = CoreState::dynamics_f(
                &state,
                &self.rac_rand_state,
                &interactions,
                world_parameters,
                parameters,
            );
            state = state + dt * delta;
        }
        let DepStates {
            geom_state,
            chem_state,
            mech_state,
        } = CoreState::calc_dep_states(&state, &self.rac_rand_state, interactions, parameters);
        state.vertex_coords = enforce_volume_exclusion(
            &state.vertex_coords,
            &geom_state,
            close_cells,
            cell_bboxes,
            cell_vcs,
        );
        // println!("++++++++++++");
        #[cfg(debug_assertions)]
        state.validate("euler", &parameters);
        // println!("{}", state);
        // let dep_vars = CellState::calc_dep_vars(&state, &self.rac_randomization, interactions, parameters);
        // println!("{}", dep_vars);
        // println!("++++++++++++");
        let rac_rand_state = match (tstep == self.rac_rand_state.next_update, rng) {
            (true, Some(cr)) => self.rac_rand_state.update(cr, tstep, parameters),
            _ => self.rac_rand_state,
        };
        //println!("{}", rac_rand_state);
        //println!("++++++++++++");
        ModelCell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_rand_state,
            geom_state,
            chem_state,
            mech_state,
        }
    }

    #[allow(unused)]
    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        interactions: &InteractionState,
        close_cells: &CloseCellInfo,
        cell_bboxes: &[Bbox],
        cell_vcs: &[&[P2D; NVERTS]],
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> ModelCell {
        // println!("using rkdp5...");
        let aux_args = AuxArgs {
            max_iters: 100,
            atol: 1e-8,
            rtol: 1e-3,
            init_h_factor: Some(0.1),
        };
        let result = rkdp5::integrator(
            1.0,
            CoreState::dynamics_f,
            &self.state,
            &self.rac_rand_state,
            interactions,
            world_parameters,
            parameters,
            aux_args,
        );

        // println!(
        //     "num_iters: {}, num_rejections: {}",
        //     result.num_iters, result.num_rejections
        // );
        let mut state = result.y.expect("too many iterations!");
        let DepStates {
            geom_state,
            chem_state,
            mech_state,
        } = CoreState::calc_dep_states(&state, &self.rac_rand_state, interactions, parameters);
        state.vertex_coords = enforce_volume_exclusion(
            &state.vertex_coords,
            &geom_state,
            close_cells,
            cell_bboxes,
            cell_vcs,
        );
        #[cfg(debug_assertions)]
        state.validate("rkdp5", parameters);
        //println!("{}", state);
        //let dep_vars = CoreState::calc_dep_vars(&state, &self.rac_rand_state, interactions, parameters);
        //println!("{}", dep_vars);
        let rac_rand_state = match (tstep == self.rac_rand_state.next_update, rng) {
            (true, Some(cr)) => self.rac_rand_state.update(cr, tstep, parameters),
            _ => self.rac_rand_state,
        };
        ModelCell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_rand_state,
            geom_state,
            chem_state,
            mech_state,
        }
    }
}

/// Calculate the area of an "ideal" initial cell of radius R, if it has n vertices.
pub fn calc_init_cell_area(r: f32, n: usize) -> f32 {
    let poly_coords = (0..n)
        .map(|vix| {
            let theta = (vix as f32) / (n as f32) * 2.0 * PI;
            P2D {
                x: r * theta.cos(),
                y: r * theta.sin(),
            }
        })
        .collect::<Vec<P2D>>();
    calc_poly_area(&poly_coords)
}
