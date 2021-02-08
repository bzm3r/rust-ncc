// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod mechanics;
pub mod py_compare;
pub mod rkdp5;
pub mod states;

use crate::cell::chemistry::RacRandState;
use crate::cell::rkdp5::AuxArgs;
use crate::cell::states::Core;
use crate::interactions::{ContactData, Interactions};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::pcg32::Pcg32;
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
        &mut self,
        tstep: u32,
        int_steps: u32,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
    ) -> Result<Cell, String> {
        let mut state = self.core;
        // Assumed normalized time by time provided in CharQuant.
        // Therefore, we can take the time period to integrate over
        // as 1.0.
        let dt = 1.0 / (int_steps as f64);
        for _ in 0..int_steps {
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
            state.enforce_volume_exclusion(
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
                tstep + 1,
                rng,
                parameters,
            ),
        })
    }

    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
    ) -> Result<Cell, String> {
        let aux_args = AuxArgs {
            max_iters: 20,
            atol: 1e-3,
            rtol: 1e-3,
            init_h_factor: Some(0.1),
        };
        let result = rkdp5::integrator(
            1.0,
            Core::derivative,
            &self.core,
            &self.rac_rand,
            interactions,
            world_parameters,
            parameters,
            aux_args,
        );

        let mut state =
            result.y.expect("rkdp5 integrator: too many iterations!");
        state
            .enforce_volume_exclusion(&self.core.poly, &contact_data)
            .map_err(|e| format!("ci={}\n{}", self.ix, e))?;

        #[cfg(feature = "validate")]
        state.validate("rkdp5")?;

        Ok(Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(tstep, rng, parameters),
        })
    }
}
