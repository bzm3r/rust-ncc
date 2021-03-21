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
pub mod states;

use crate::cell::chemistry::RacRandState;
use crate::cell::states::{confirm_volume_exclusion, Core};
use crate::interactions::{ContactData, Interactions};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::pcg32::Pcg32;
use crate::world::{EulerOpts, RkOpts};
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
    /// State of Random Rac1 activity that affected `core`.
    pub rac_rand: RacRandState,
    /// Core state of the cell (position, Rho GTPase).
    pub core: Core,
}

impl Cell {
    pub fn new(
        ix: usize,
        group_ix: usize,
        core: Core,
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

    /// Suppose our current state is `state`. We want to determine
    /// the next state after a time period `dt` has elapsed. We
    /// assume `(next_state - state)/delta(t) = delta(state)`.
    pub fn simulate_euler(
        &self,
        tpoint: f64,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
        int_opts: EulerOpts,
    ) -> Result<Cell, String> {
        let mut state = self.core;
        // Assumed normalized time by time provided in CharQuant.
        // Therefore, we can take the time period to integrate over
        // as 1.0.
        let dt = 1.0 / (int_opts.num_int_steps as f64);
        for _ in 0..int_opts.num_int_steps {
            // d(state)/dt = dynamics_f(state) <- calculate RHS of ODE
            let delta = Core::derivative(
                &state,
                &self.rac_rand,
                &interactions,
                world_parameters,
                parameters,
            );
            state = state + delta.time_step(dt);
            // Enforcing volume exclusion! Tricky!
            state.strict_enforce_volume_exclusion(
                self.ix,
                &self.core.poly,
                &contact_data,
            )?;
        }

        #[cfg(feature = "validate")]
        state.validate("euler")?;

        Ok(Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(
                tpoint + 1.0,
                rng,
                parameters,
            ),
        })
    }

    /// Suppose our current state is `state`. We want to determine
    /// the next state after a time period `dt` has elapsed. We
    /// assume `(next_state - state)/delta(t) = delta(state)`.
    pub fn simulate_euler_debug(
        &self,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        int_opts: EulerOpts,
    ) -> Result<Vec<Cell>, String> {
        let cell_ix = self.ix;
        // println!("cell_ix: {}", cell_ix);
        let mut r: Vec<Cell> =
            Vec::with_capacity(int_opts.num_int_steps as usize);
        let mut state = self.core;
        let dt = 1.0 / (int_opts.num_int_steps as f64);
        confirm_volume_exclusion(
            cell_ix,
            &self.core.poly,
            &contact_data,
            "old_vs",
        )?;
        // let (focus_vi, other_focus_v) = if cell_ix == 0 {
        //     (0, contact_data[0].poly.verts[8])
        // } else {
        //     (8, contact_data[0].poly.verts[0])
        // };
        for _ in 0..int_opts.num_int_steps {
            let delta = Core::derivative(
                &state,
                &self.rac_rand,
                &interactions,
                world_parameters,
                parameters,
            );
            let old_vs = state.poly;
            state = state + delta.time_step(dt);
            // println!(
            //     "before vol_ex: state.poly[{}] = {} | other: {}",
            //     focus_vi, state.poly[focus_vi], other_focus_v
            // );
            // Enforcing volume exclusion! Tricky!
            state.strict_enforce_volume_exclusion(
                self.ix,
                &old_vs,
                &contact_data,
            )?;
            // println!(
            //     "after vol_ex | state.poly[{}] = {} | other: {}",
            //     focus_vi, state.poly[focus_vi], other_focus_v
            // );
            r.push(Cell {
                ix: self.ix,
                group_ix: self.group_ix,
                core: state,
                rac_rand: self.rac_rand,
            })
        }
        #[cfg(feature = "validate")]
        state.validate("euler")?;
        // println!(
        //     "final | state.poly[{}] = {} | other: {}",
        //     focus_vi,
        //     r.last().unwrap().core.poly[focus_vi],
        //     other_focus_v
        // );
        // println!(
        //     "final poly: {}",
        //     poly_to_string(&r.last().unwrap().core.poly)
        // );
        // println!(
        //     "other poly: {}",
        //     poly_to_string(&contact_data[0].poly.verts)
        // );
        Ok(r)
    }

    pub fn simulate_rkdp5(
        &self,
        tpoint: f64,
        dt: f64,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
        int_opts: RkOpts,
    ) -> Result<Cell, String> {
        let mut result = rkdp5::integrator(
            dt,
            Core::derivative,
            self.core,
            &self.rac_rand,
            interactions,
            world_parameters,
            parameters,
            &contact_data,
            int_opts,
        );

        match &mut result.state {
            Ok(cs) => {
                #[cfg(feature = "validate")]
                cs.validate("rkdp5")?;

                Ok(Cell {
                    ix: self.ix,
                    group_ix: self.group_ix,
                    core: *cs,
                    rac_rand: self
                        .rac_rand
                        .update(tpoint, rng, parameters),
                })
            }
            Err(err) => Err(format!("{} ({:?})", err, int_opts)),
        }
    }
}
