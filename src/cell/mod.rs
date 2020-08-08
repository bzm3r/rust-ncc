// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod mechanics;
pub mod rkdp5;
pub mod state;

use crate::cell::chemistry::{gen_rgtp_distrib, RacRandState, RgtpLayout};
use crate::cell::rkdp5::AuxArgs;
use crate::cell::state::State;
use crate::math::geometry::calc_poly_area;
use crate::math::p2d::P2D;
use crate::parameters::Parameters;
use crate::world::interactions::InteractionState;
use crate::world::RandomEventGenerator;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

/// Model cell.
#[derive(Copy, Clone, Deserialize, Serialize, Schematize)]
pub struct Cell {
    /// Cell index.
    pub ix: u32,
    /// Index of cell type.
    pub group_ix: u32,
    pub state: State,
    pub rac_rand_state: RacRandState,
}

fn gen_vertex_coords(centroid: P2D, radius: f32) -> [P2D; NVERTS as usize] {
    let mut r = [P2D::default(); NVERTS as usize];
    (0..NVERTS as usize).for_each(|vix| {
        let vf = (vix as f32) / (NVERTS as f32);
        let theta = 2.0 * PI * vf;
        r[vix] = P2D {
            x: centroid.x + theta.cos() * radius,
            y: centroid.y + theta.sin() * radius,
        };
    });
    r
}

#[allow(unused)]
pub enum VertexGenInfo {
    Centroid(P2D),
    Coordinates([P2D; NVERTS as usize]),
}

impl Cell {
    pub fn new(
        ix: u32,
        group_ix: u32,
        vertex_gen_info: VertexGenInfo,
        parameters: &Parameters,
        reg: Option<&mut RandomEventGenerator>,
    ) -> Cell {
        let vertex_coords = match vertex_gen_info {
            VertexGenInfo::Centroid(centroid) => gen_vertex_coords(centroid, parameters.cell_r),
            VertexGenInfo::Coordinates(vcs) => vcs,
        };
        let (rac_acts, rac_inacts) = gen_rgtp_distrib(
            parameters.init_frac_active,
            parameters.init_frac_inactive,
            &RgtpLayout::Random,
            //&RgtpLayout::BiasedVertices(vec![0, 1, 2, 3]),
        );
        let (rho_acts, rho_inacts) = gen_rgtp_distrib(
            parameters.init_frac_active,
            parameters.init_frac_inactive,
            &RgtpLayout::Random,
            //&RgtpLayout::BiasedVertices(vec![8, 9, 10, 11]),
        );
        let state = State::new(vertex_coords, rac_acts, rac_inacts, rho_acts, rho_inacts);
        let rac_rand_state = RacRandState::init(
            match reg {
                Some(cr) => Some(&mut cr.rng),
                None => None,
            },
            parameters,
        );
        Cell {
            ix,
            group_ix,
            state,
            rac_rand_state,
        }
    }

    #[allow(unused)]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        inter_state: &InteractionState,
        rng: Option<&mut RandomEventGenerator>,
        parameters: &Parameters,
    ) -> Cell {
        let mut state = self.state;
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for i in 0..nsteps {
            println!("++++++++++++");
            println!("{}", state);
            let dep_vars =
               State::calc_dep_vars(&state, &self.rac_rand_state, inter_state, parameters);
            println!("{}", dep_vars);
            println!("++++++++++++");
            let delta = State::dynamics_f(&state, &self.rac_rand_state, &inter_state, parameters);
            state = state + dt * delta;
        }
        // println!("++++++++++++");
        #[cfg(debug_assertions)]
        state.validate("euler", &parameters);
        // println!("{}", state);
        // let dep_vars = CellState::calc_dep_vars(&state, &self.rac_randomization, inter_state, parameters);
        // println!("{}", dep_vars);
        // println!("++++++++++++");
        let rac_rand_state = match (tstep == self.rac_rand_state.next_update, rng) {
            (true, Some(cr)) => self.rac_rand_state.update(cr, tstep, parameters),
            _ => self.rac_rand_state,
        };
        println!("{}", rac_rand_state);
        println!("++++++++++++");
        Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_rand_state,
        }
    }

    #[allow(unused)]
    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        inter_state: &InteractionState,
        rng: Option<&mut RandomEventGenerator>,
        parameters: &Parameters,
    ) -> Cell {
        // println!("using rkdp5...");
        let aux_args = AuxArgs {
            max_iters: 100,
            atol: 1e-8,
            rtol: 1e-3,
            init_h_factor: Some(0.1),
        };
        let result = rkdp5::integrator(
            1.0,
            State::dynamics_f,
            &self.state,
            &self.rac_rand_state,
            inter_state,
            parameters,
            aux_args,
        );

        // println!(
        //     "num_iters: {}, num_rejections: {}",
        //     result.num_iters, result.num_rejections
        // );
        let state = result.y.expect("too many iterations!");
        #[cfg(debug_assertions)]
        state.validate("rkdp5", parameters);
        //println!("{}", state);
        //let dep_vars = State::calc_dep_vars(&state, &self.rac_rand_state, inter_state, parameters);
        //println!("{}", dep_vars);
        let rac_rand_state = match (tstep == self.rac_rand_state.next_update, rng) {
            (true, Some(cr)) => self.rac_rand_state.update(cr, tstep, parameters),
            _ => self.rac_rand_state,
        };
        Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_rand_state,
        }
    }
}

/// Calculate the area of an "ideal" initial cell of radius R, if it has n vertices.
pub fn calc_init_cell_area(r: f32, n: u32) -> f32 {
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
