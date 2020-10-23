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

use crate::interactions::CellInteractions;
use crate::math::geometry::{calc_poly_area, is_point_in_poly, BBox};
use crate::math::p2d::V2D;
use crate::model_cell::chemistry::RacRandState;
use crate::model_cell::core_state::{ChemState, CoreState, GeomState, MechState};
use crate::model_cell::rkdp5::AuxArgs;
use crate::parameters::{GlobalParameters, Parameters};
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

fn move_point_out(
    mut out_p: V2D,
    mut in_p: V2D,
    poly_bbox: &BBox,
    poly: &[V2D],
    num_iters: usize,
) -> V2D {
    let mut n = 0;
    while n < num_iters {
        let new_p = 0.5 * (out_p + in_p);
        if is_point_in_poly(&new_p, poly_bbox, poly) {
            in_p = new_p;
        } else {
            out_p = new_p;
        }
        n += 1;
    }
    out_p
}

#[cfg(debug_assertions)]
pub fn confirm_volume_exclusion(
    vcs: &[V2D; NVERTS],
    contact_polys: &[(BBox, [V2D; NVERTS])],
    panic_label: &str,
) {
    for (vi, vc) in vcs.iter().enumerate() {
        for (bbox, poly) in contact_polys {
            if is_point_in_poly(vc, bbox, poly) {
                panic!(
                    "{:?} violate volume exclusion. vcs[{:?}] = {:?}, other poly = {:?}",
                    panic_label, vi, vc, poly
                );
            }
        }
    }
}

fn enforce_volume_exclusion(
    old_vcs: &[V2D; NVERTS],
    mut new_vcs: [V2D; NVERTS],
    contact_polys: Vec<(BBox, [V2D; NVERTS])>,
) -> [V2D; NVERTS] {
    #[cfg(debug_assertions)]
    confirm_volume_exclusion(old_vcs, contact_polys.as_slice(), "old_vcs");
    for (old_vc, new_vc) in old_vcs.iter().zip(new_vcs.iter_mut()) {
        for (obb, ovcs) in contact_polys.iter() {
            if is_point_in_poly(new_vc, obb, ovcs) {
                let out_p = *old_vc;
                let in_p = *new_vc;
                *new_vc = move_point_out(out_p, in_p, obb, ovcs, 3);
            }
        }
    }
    #[cfg(debug_assertions)]
    confirm_volume_exclusion(&new_vcs, &contact_polys, "new_vcs");

    new_vcs
}

impl ModelCell {
    pub fn new(
        ix: u32,
        group_ix: u32,
        state: CoreState,
        interactions: &CellInteractions,
        parameters: &Parameters,
        reg: Option<&mut RandomEventGenerator>,
    ) -> ModelCell {
        let rac_rand_state = RacRandState::init(
            match reg {
                Some(cr) => Some(&mut cr.rng),
                None => None,
            },
            parameters,
        );
        let geom_state = state.calc_geom_state();
        let mech_state = state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &rac_rand_state,
            &interactions,
            parameters,
        );
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
        interactions: &CellInteractions,
        contact_polys: Vec<(BBox, [V2D; NVERTS])>,
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &GlobalParameters,
        parameters: &Parameters,
    ) -> ModelCell {
        let mut state = self.state;
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for i in 0..nsteps {
            //println!("++++++++++++");
            //println!("{}", state);
            //let dep_vars =
            //    CoreState::calc_dep_states(&state, &self.rac_rand_state, interactions, parameters);
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
        let geom_state = state.calc_geom_state();
        let mech_state = state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand_state,
            &interactions,
            parameters,
        );
        state.vertex_coords = enforce_volume_exclusion(
            &self.state.vertex_coords,
            state.vertex_coords,
            contact_polys,
        );
        let geom_state = state.calc_geom_state();
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
        interactions: &CellInteractions,
        contact_polys: Vec<(BBox, [V2D; NVERTS])>,
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &GlobalParameters,
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
        let geom_state = state.calc_geom_state();
        let mech_state = state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand_state,
            &interactions,
            parameters,
        );
        let mv_dirs = mech_state
            .sum_fs
            .iter()
            .map(|sf| -1.0 * sf.unitize())
            .collect::<Vec<V2D>>();
        state.vertex_coords = enforce_volume_exclusion(
            &self.state.vertex_coords,
            state.vertex_coords,
            contact_polys,
        );
        let geom_state = state.calc_geom_state();
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
            V2D {
                x: r * theta.cos(),
                y: r * theta.sin(),
            }
        })
        .collect::<Vec<V2D>>();
    calc_poly_area(&poly_coords)
}
