// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod core_state;
// pub mod geometry;
pub mod mechanics;
pub mod rkdp5;

use crate::cell::chemistry::RacRandState;
use crate::cell::core_state::{
    ChemState, CoreState, GeomState, MechState,
};
use crate::cell::rkdp5::AuxArgs;
use crate::interactions::{CellInteractions, ContactData};
use crate::math::geometry::{calc_poly_area, LineSeg2D};
use crate::math::v2d::{poly_to_string, V2D};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

/// Cell state structure.
#[derive(Copy, Clone, Deserialize, Serialize, Schematize)]
pub struct CellState {
    /// Index of cell within world.
    pub ix: u32,
    /// Index of group that cell belongs to.
    pub group_ix: u32,
    pub core: CoreState,
    pub rac_rand: RacRandState,
    pub mech: MechState,
    pub geom: GeomState,
    pub chem: ChemState,
}

pub fn violates_volume_exclusion(
    test_v: &V2D,
    test_u: &V2D,
    test_w: &V2D,
    contacts: &[ContactData],
) -> bool {
    let lsegs = [
        LineSeg2D::new(test_u, test_v),
        LineSeg2D::new(test_v, test_w),
    ];
    for lseg in lsegs.iter() {
        for cd in contacts.iter() {
            if lseg.intersects_poly(&cd.poly) {
                return true;
            }
        }
    }
    false
}

fn move_point_out(
    new_u: &V2D,
    mut new_v: V2D,
    new_w: &V2D,
    mut good_v: V2D,
    contacts: &[ContactData],
    num_iters: usize,
) -> V2D {
    let mut n = 0;
    while n < num_iters {
        let test_v = 0.5 * (good_v + new_v);
        if violates_volume_exclusion(&test_v, new_u, new_w, contacts)
        {
            new_v = test_v;
        } else {
            good_v = test_v;
        }
        n += 1;
    }
    good_v
}

#[cfg(feature = "debug_mode")]
pub fn confirm_volume_exclusion(
    vs: &[V2D; NVERTS],
    contacts: &[ContactData],
    msg: &str,
) -> Result<(), String> {
    for (vi, v) in vs.iter().enumerate() {
        let u = &vs[circ_ix_minus(vi, NVERTS)];
        let w = &vs[circ_ix_plus(vi, NVERTS)];
        for ContactData { poly, .. } in contacts {
            if violates_volume_exclusion(v, u, w, contacts) {
                return Err(format!(
                    "{} violates volume exclusion.\n\
                    vs[{}] = {}, \n\
                    other poly = {}",
                    msg,
                    vi,
                    v,
                    &poly_to_string(&poly.verts)
                ));
            }
        }
    }
    Ok(())
}

#[cfg(feature = "debug_mode")]
pub fn enforce_volume_exclusion(
    old_vs: &[V2D; NVERTS],
    mut new_vs: [V2D; NVERTS],
    contacts: Vec<ContactData>,
) -> Result<[V2D; NVERTS], String> {
    confirm_volume_exclusion(&old_vs, &contacts, "old_vs")?;
    for vi in 0..NVERTS {
        let old_v = old_vs[vi];
        let new_v = new_vs[vi];
        let new_u = new_vs[circ_ix_minus(vi, NVERTS)];
        let new_w = new_vs[circ_ix_plus(vi, NVERTS)];
        new_vs[vi] = move_point_out(
            &new_u, new_v, &new_w, old_v, &contacts, 5,
        );
    }
    confirm_volume_exclusion(&new_vs, &contacts, "new_vs")?;
    Ok(new_vs)
}

#[cfg(not(feature = "debug_mode"))]
pub fn enforce_volume_exclusion(
    old_vs: &[V2D; NVERTS],
    mut new_vs: [V2D; NVERTS],
    contact_polys: Vec<(BBox, [V2D; NVERTS])>,
) -> [V2D; NVERTS] {
    for (old_v, new_v) in old_vs.iter().zip(new_vs.iter_mut()) {
        for (obb, ovs) in contact_polys.iter() {
            if is_point_in_poly(new_v, obb, ovs) {
                *new_v = move_point_out(*old_v, *new_v, obb, ovs, 3);
            }
        }
    }
    new_vs
}

impl CellState {
    pub fn new(
        ix: u32,
        group_ix: u32,
        state: CoreState,
        interactions: &CellInteractions,
        parameters: &Parameters,
        rac_rand_state: RacRandState,
    ) -> CellState {
        let geom_state = state.calc_geom_state();
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &rac_rand_state,
            &interactions,
            parameters,
        );
        CellState {
            ix,
            group_ix,
            core: state,
            rac_rand: rac_rand_state,
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        }
    }

    #[allow(unused)]
    #[cfg(feature = "debug_mode")]
    /// Suppose our current state is `state`. We want to determine
    /// the next state after a time period `dt` has elapsed. We
    /// assume `(next_state - state)/delta(t) = delta(state)`.
    pub fn simulate_euler(
        &mut self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> Result<CellState, String> {
        let mut state = self.core;
        let nsteps: u32 = 10;
        // Assumed normalized time by time provided in CharQuant.
        // Therefore, we can take the time period to integrate over
        // as 1.0.
        let dt = 1.0 / (nsteps as f32);
        for _ in 0..nsteps {
            // d(state)/dt = dynamics_f(state) <- calculate RHS of ODE
            let delta = CoreState::dynamics_f(
                &state,
                &self.rac_rand,
                &interactions,
                world_parameters,
                parameters,
            );
            state = state + dt * delta;
        }
        let geom_state = state.calc_geom_state();
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        // Enforcing volume exclusion! Tricky!
        state.poly = enforce_volume_exclusion(
            &self.core.poly,
            state.poly,
            contact_data,
        )?;
        let geom_state = state.calc_geom_state();
        // println!("++++++++++++");
        state.validate("euler", &parameters)?;
        Ok(CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(tstep, parameters),
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        })
    }

    #[allow(unused)]
    #[cfg(not(feature = "debug_mode"))]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_polys: Vec<(BBox, [V2D; NVERTS])>,
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> CellState {
        let mut state = self.core;
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for i in 0..nsteps {
            let delta = CoreState::dynamics_f(
                &state,
                &self.rac_rand,
                &interactions,
                world_parameters,
                parameters,
            );
            state = state + dt * delta;
        }
        let geom_state = state.calc_geom_state();
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        state.poly = enforce_volume_exclusion(
            &self.core.poly,
            state.poly,
            contact_polys,
        );
        let geom_state = state.calc_geom_state();
        // println!("++++++++++++");
        let rac_rand_state =
            match (tstep == self.rac_rand.next_update, rng) {
                (true, Some(cr)) => {
                    self.rac_rand.update(cr, tstep, parameters)
                }
                _ => self.rac_rand,
            };
        CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: rac_rand_state,
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        }
    }

    #[cfg(feature = "debug_mode")]
    pub fn simulate_rkdp5(
        &mut self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> Result<CellState, String> {
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
            &self.core,
            &self.rac_rand,
            interactions,
            world_parameters,
            parameters,
            aux_args,
        );

        // println!(
        //     "num_iters: {}, num_rejections: {}",
        //     result.num_iters, result.num_rejections
        // );
        let mut state =
            result.y.expect("rkdp5 integrator: too many iterations!");
        let geom_state = state.calc_geom_state();
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        state.poly = enforce_volume_exclusion(
            &self.core.poly,
            state.poly,
            contact_data,
        )
        .map_err(|e| format!("ci={}\n{}", self.ix, e))?;
        let geom_state = state.calc_geom_state();
        state.validate("rkdp5", &parameters)?;
        // println!("{}", state);
        //let dep_vars = CoreState::calc_dep_vars(&state, &self.rac_rand_state, interactions, parameters);
        // println!("{}", dep_vars);
        Ok(CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(tstep, parameters),
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        })
    }

    #[cfg(not(feature = "debug_mode"))]
    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_polys: Vec<(BBox, [V2D; NVERTS])>,
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> CellState {
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
            &self.core,
            &self.rac_rand,
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
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
            &geom_state,
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        let mv_dirs = mech_state
            .sum_fs
            .iter()
            .map(|sf| -1.0 * sf.unitize())
            .collect::<Vec<V2D>>();
        state.poly = enforce_volume_exclusion(
            &self.core.poly,
            state.poly,
            contact_polys,
        );
        let geom_state = state.calc_geom_state();
        #[cfg(feature = "debug_mode")]
        state.validate("rkdp5", parameters);
        // println!("{}", state);
        //let dep_vars = CoreState::calc_dep_vars(&state, &self.rac_rand_state, interactions, parameters);
        // println!("{}", dep_vars);
        let rac_rand_state =
            match (tstep == self.rac_rand.next_update, rng) {
                (true, Some(cr)) => {
                    self.rac_rand.update(cr, tstep, parameters)
                }
                _ => self.rac_rand,
            };
        CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: rac_rand_state,
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
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
