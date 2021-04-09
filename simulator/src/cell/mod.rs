// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod delta_cell;
pub mod mechanics;
pub mod rkdp5;

use crate::cell::chemistry::RacRandState;
use crate::cell::delta_cell::{
    confirm_volume_exclusion, Core, GeomState,
};
use crate::interactions::{
    Contact, InteractionGenerator, Interactions,
};
use crate::math::v2d::V2d;
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::pcg32::Pcg32;
use crate::world::{EulerOpts, RkOpts};
use crate::NVERTS;
use serde::{Deserialize, Serialize};

/// Cell state structure.
#[derive(
    Copy, Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct Cell {
    /// Index of cell within world.
    pub ix: usize,
    /// Index of group that cell belongs to.
    pub group_ix: usize,
    /// State of Random Rac1 activity that affected this cell.
    pub rac_rand: RacRandState,
    /// Polygon representing cell shape.
    pub poly: [V2d; NVERTS],
    /// Fraction of Rac1 active at each vertex.
    pub rac_acts: [f64; NVERTS],
    /// Fraction of Rac1 inactive at each vertex.
    pub rac_inacts: [f64; NVERTS],
    /// Fraction of RhoA active at each vertex.
    pub rho_acts: [f64; NVERTS],
    /// Fraction of RhoA inactive at each vertex.
    pub rho_inacts: [f64; NVERTS],
    /// COA factor at each vertex.
    pub x_coas: [f64; NVERTS],
    /// CIL factor at each vertex.
    pub x_cils: [f64; NVERTS],
    /// CAL factor at each vertex,
    pub x_cals: [f64; NVERTS],
    /// Chemoattractant factor at each vertex,
    pub x_chemoas: [f64; NVERTS],
    /// Adhesive force acting at each vertex,
    pub x_adhs: [V2d; NVERTS],
    /// Geometric state due to vertex positions.
    pub geom: GeomState,
}

impl Cell {
    pub fn new(
        ix: usize,
        group_ix: usize,
        parameters: &Parameters,
        rng: &mut Pcg32,
    ) -> Cell {
        let rac_rand = if parameters.randomization {
            RacRandState::new(rng, parameters)
        } else {
            RacRandState::default()
        };
        Cell {
            ix,
            group_ix,
            core,
            rac_rand,
        }
    }

    pub fn update_core(
        &self,
        poly: [V2d; NVERTS],
        rac_acts: [f64; NVERTS],
        rac_inacts: [f64; NVERTS],
        rho_acts: [f64; NVERTS],
        rho_inacts: [f64; NVERTS],
    ) -> Cell {
        Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            rac_rand: self.rac_rand,
            poly,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
            x_coas: self.x_coas,
            x_cils: self.x_cils,
            x_cals: self.x_cals,
            x_chemoas: self.x_chemoas,
            x_adhs: self.x_adhs,
            geom: self.geom,
        }
    }

    // /// Suppose our current state is `state`. We want to determine
    // /// the next state after a time period `dt` has elapsed. We
    // /// assume `(next_state - state)/delta(t) = delta(state)`.
    // pub fn simulate_euler(
    //     &self,
    //     tpoint: f64,
    //     interactions: &Interactions,
    //     contact_data: Vec<Contact>,
    //     world_parameters: &WorldParameters,
    //     parameters: &Parameters,
    //     rng: &mut Pcg32,
    //     int_opts: EulerOpts,
    // ) -> Result<Cell, String> {
    //     let mut state = self.core;
    //     // Assumed normalized time by time provided in CharQuant.
    //     // Therefore, we can take the time period to integrate over
    //     // as 1.0.
    //     let dt = 1.0 / (int_opts.num_int_steps as f64);
    //     for _ in 0..int_opts.num_int_steps {
    //         // d(state)/dt = dynamics_f(state) <- calculate RHS of ODE
    //         let delta = Core::derivative(
    //             &state,
    //             &self.rac_rand,
    //             &interactions,
    //             world_parameters,
    //             parameters,
    //             x_adhs,
    //             x_coas,
    //             x_cils,
    //             x_cals,
    //             x_chemoas,
    //         );
    //         state = state + delta.time_step(dt);
    //         // Enforcing volume exclusion! Tricky!
    //         state.strict_enforce_volume_exclusion(
    //             &self.core.poly,
    //             &contact_data,
    //         )?;
    //     }
    //
    //     #[cfg(feature = "validate")]
    //     state.validate("euler")?;
    //
    //     Ok(Cell {
    //         ix: self.ix,
    //         group_ix: self.group_ix,
    //         core: state,
    //         rac_rand: self.rac_rand.update(
    //             tpoint + 1.0,
    //             rng,
    //             parameters,
    //         ),
    //     })
    // }

    // /// Suppose our current state is `state`. We want to determine
    // /// the next state after a time period `dt` has elapsed. We
    // /// assume `(next_state - state)/delta(t) = delta(state)`.
    // pub fn simulate_euler_debug(
    //     &self,
    //     interactions: &Interactions,
    //     contact_data: Vec<Contact>,
    //     world_parameters: &WorldParameters,
    //     parameters: &Parameters,
    //     int_opts: EulerOpts,
    // ) -> Result<Vec<Cell>, String> {
    //     // println!("cell_ix: {}", cell_ix);
    //     let mut r: Vec<Cell> =
    //         Vec::with_capacity(int_opts.num_int_steps as usize);
    //     let mut state = self.core;
    //     let dt = 1.0 / (int_opts.num_int_steps as f64);
    //     confirm_volume_exclusion(
    //         &self.core.poly,
    //         &contact_data,
    //         "old_vs",
    //     )?;
    //     // let (focus_vi, other_focus_v) = if cell_ix == 0 {
    //     //     (0, contact_data[0].poly.verts[8])
    //     // } else {
    //     //     (8, contact_data[0].poly.verts[0])
    //     // };
    //     for _ in 0..int_opts.num_int_steps {
    //         let delta = Core::derivative(
    //             &state,
    //             &self.rac_rand,
    //             &interactions,
    //             world_parameters,
    //             parameters,
    //         );
    //         let old_vs = state.poly;
    //         state = state + delta.time_step(dt);
    //         // println!(
    //         //     "before vol_ex: state.poly[{}] = {} | other: {}",
    //         //     focus_vi, state.poly[focus_vi], other_focus_v
    //         // );
    //         // Enforcing volume exclusion! Tricky!
    //         state.strict_enforce_volume_exclusion(
    //             &old_vs,
    //             &contact_data,
    //         )?;
    //         // println!(
    //         //     "after vol_ex | state.poly[{}] = {} | other: {}",
    //         //     focus_vi, state.poly[focus_vi], other_focus_v
    //         // );
    //         r.push(Cell {
    //             ix: self.ix,
    //             group_ix: self.group_ix,
    //             core: state,
    //             rac_rand: self.rac_rand,
    //         })
    //     }
    //     #[cfg(feature = "validate")]
    //     state.validate("euler")?;
    //     // println!(
    //     //     "final | state.poly[{}] = {} | other: {}",
    //     //     focus_vi,
    //     //     r.last().unwrap().core.poly[focus_vi],
    //     //     other_focus_v
    //     // );
    //     // println!(
    //     //     "final poly: {}",
    //     //     poly_to_string(&r.last().unwrap().core.poly)
    //     // );
    //     // println!(
    //     //     "other poly: {}",
    //     //     poly_to_string(&contact_data[0].poly.verts)
    //     // );
    //     Ok(r)
    // }

    pub fn simulate_rkdp5(
        &self,
        tpoint: f64,
        dt: f64,
        inter_gen: &mut InteractionGenerator,
        contacts: Vec<Contact>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
        int_opts: RkOpts,
    ) -> Result<Cell, String> {
        let mut result = rkdp5::integrate(
            dt,
            Core::derivative,
            self.core,
            &self.rac_rand,
            inter_gen,
            world_parameters,
            parameters,
            int_opts,
        );

        match &mut result.state {
            Ok(cs) => {
                #[cfg(feature = "validate")]
                cs.validate("rkdp5")?;

                cs.strict_enforce_volume_exclusion(
                    &self.core.poly,
                    &contacts,
                )?;

                Ok(Cell {
                    ix: self.ix,
                    group_ix: self.group_ix,
                    core: *cs,
                    rac_rand: self
                        .rac_rand
                        .update(tpoint, rng, parameters),
                })
            }
            Err(e) => Err(e.into()),
        }
    }
}
