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
use crate::math::geometry::{
    calc_poly_area, is_point_in_poly, BBox, LineSeg2D,
};
use crate::math::v2d::{poly_to_string, V2D};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::world::RandomEventGenerator;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

/// Cell state structure.
#[derive(Copy, Clone, Deserialize, Serialize, Schematize)]
pub struct CellState {
    /// Cell index.
    pub ix: u32,
    /// Index of cell type.
    pub group_ix: u32,
    pub core: CoreState,
    pub rac_rand: RacRandState,
    pub mech: MechState,
    pub geom: GeomState,
    pub chem: ChemState,
}

pub fn check_poly_intersect(
    u: &V2D,
    v: &V2D,
    w: &V2D,
    poly_bbox: &BBox,
    poly: &[V2D; NVERTS],
) -> bool {
    let ls_uv = LineSeg2D::new(u, v);
    let ls_vw = LineSeg2D::new(v, w);
    let opt_bb = Some(poly_bbox);
    let uv_intersects = ls_uv.intersects_poly(opt_bb, poly);
    let vw_intersects = ls_vw.intersects_poly(opt_bb, poly);
    let is_v_in_poly = is_point_in_poly(v, Some(poly_bbox), poly);

    uv_intersects || vw_intersects || is_v_in_poly
}

fn move_point_out(
    bad_vi: usize,
    good_vs: &[V2D; NVERTS],
    mut bad_vs: [V2D; NVERTS],
    // new_u: V2D,
    // new_w: V2D,
    close_verts: &[usize],
    poly_bbox: &BBox,
    poly: &[V2D; NVERTS],
    num_iters: usize,
) -> V2D {
    let mut good_v = good_vs[bad_vi];
    let mut n = 0;
    while n < num_iters {
        let test_v = 0.5 * (good_v + bad_vs[bad_vi]);
        // check_poly_intersect(
        //     &new_u, &test_v, &new_w, poly_bbox, poly,
        // )
        if is_point_in_poly(&test_v, Some(poly_bbox), poly)
            || close_verts.iter().any(|&ovi| {
                is_point_in_poly(
                    &poly[ovi],
                    Some(&BBox::from_points(&bad_vs)),
                    &bad_vs,
                )
            })
        {
            bad_vs[bad_vi] = test_v;
        } else {
            good_v = test_v;
        }
        n += 1;
    }
    good_v
}

#[cfg(feature = "custom_debug")]
pub fn confirm_volume_exclusion(
    vs: &[V2D; NVERTS],
    contact_data: &[ContactData],
    msg: &str,
) -> Result<(), String> {
    for (vi, v) in vs.iter().enumerate() {
        for ContactData { poly, poly_bb, .. } in contact_data {
            if is_point_in_poly(v, Some(poly_bb), poly) {
                return Err(format!(
                    "{} violates volume exclusion.\n\
                    vs[{}] = {}, \n\
                    other poly = {}",
                    msg,
                    vi,
                    v,
                    &poly_to_string(poly)
                ));
            }
        }
    }
    Ok(())
}

#[cfg(feature = "custom_debug")]
pub fn enforce_volume_exclusion(
    old_vs: &[V2D; NVERTS],
    mut new_vs: [V2D; NVERTS],
    contact_data: Vec<ContactData>,
) -> Result<[V2D; NVERTS], String> {
    confirm_volume_exclusion(&old_vs, &contact_data, "old_vs")?;
    for vi in 0..NVERTS {
        let nv = new_vs[vi];
        // let nw = new_vs[circ_ix_plus(vi, NVERTS)];
        // let nu = new_vs[circ_ix_minus(vi, NVERTS)];
        for ContactData {
            oci,
            close_verts,
            poly,
            poly_bb,
        } in contact_data.iter()
        {
            // check_poly_intersect(&nu, &nv, &nw, obb, opoly)
            if is_point_in_poly(&nv, Some(poly_bb), poly) {
                new_vs[vi] = move_point_out(
                    vi,
                    old_vs,
                    new_vs,
                    close_verts,
                    poly_bb,
                    poly,
                    5,
                );
            }
        }
    }
    confirm_volume_exclusion(&new_vs, &contact_data, "new_vs")?;
    Ok(new_vs)
}

#[cfg(not(feature = "custom_debug"))]
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
        reg: Option<&mut RandomEventGenerator>,
    ) -> CellState {
        let rac_rand_state = RacRandState::init(
            match reg {
                Some(cr) => Some(&mut cr.rng),
                None => None,
            },
            parameters,
        );
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
    #[cfg(feature = "custom_debug")]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_data: Vec<ContactData>,
        rng: Option<&mut RandomEventGenerator>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> Result<CellState, String> {
        let mut state = self.core;
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for _ in 0..nsteps {
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
        state.vertex_coords = enforce_volume_exclusion(
            &self.core.vertex_coords,
            state.vertex_coords,
            contact_data,
        )?;
        let geom_state = state.calc_geom_state();
        // println!("++++++++++++");
        state.validate("euler", &parameters)?;
        let rac_rand_state =
            match (tstep == self.rac_rand.next_update, rng) {
                (true, Some(cr)) => {
                    self.rac_rand.update(cr, tstep, parameters)
                }
                _ => self.rac_rand,
            };
        Ok(CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: rac_rand_state,
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        })
    }

    #[allow(unused)]
    #[cfg(not(feature = "custom_debug"))]
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
        state.vertex_coords = enforce_volume_exclusion(
            &self.core.vertex_coords,
            state.vertex_coords,
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

    #[cfg(feature = "custom_debug")]
    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        interactions: &CellInteractions,
        contact_data: Vec<ContactData>,
        rng: Option<&mut RandomEventGenerator>,
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
        state.vertex_coords = enforce_volume_exclusion(
            &self.core.vertex_coords,
            state.vertex_coords,
            contact_data,
        )
        .map_err(|e| format!("ci={}\n{}", self.ix, e))?;
        let geom_state = state.calc_geom_state();
        state.validate("rkdp5", &parameters)?;
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
        Ok(CellState {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: rac_rand_state,
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        })
    }

    #[cfg(not(feature = "custom_debug"))]
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
        state.vertex_coords = enforce_volume_exclusion(
            &self.core.vertex_coords,
            state.vertex_coords,
            contact_polys,
        );
        let geom_state = state.calc_geom_state();
        #[cfg(feature = "custom_debug")]
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
