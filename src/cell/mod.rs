// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
pub mod chemistry;
pub mod states;
// pub mod geometry;
pub mod mechanics;
pub mod pycomp;
pub mod rkdp5;

use crate::cell::chemistry::RacRandState;
use crate::cell::rkdp5::AuxArgs;
use crate::cell::states::{ChemState, Core, GeomState, MechState};
use crate::interactions::{ContactData, Interactions};
use crate::math::close_to_zero;
use crate::math::geometry::{
    calc_poly_area, check_strong_intersection,
};
use crate::math::v2d::V2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::pcg32::Pcg32;
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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
    /// Geometry (unit inward vecs, etc.) due to `core`.
    pub geom: GeomState,
    /// State of Mechanical activity (forces) due to `core``.
    pub mech: MechState,
    /// Chemical state (various reaction rates) due to `core`.
    pub chem: ChemState,
}

pub fn violates_volume_exclusion(
    test_v: &V2D,
    test_u: &V2D,
    test_w: &V2D,
    contacts: &[ContactData],
) -> bool {
    let point_pairs = [(test_u, test_v), (test_v, test_w)];
    for (p0, p1) in point_pairs.iter() {
        for cd in contacts.iter() {
            for edge in cd.poly.edges.iter() {
                if check_strong_intersection(p0, p1, edge) {
                    return true;
                }
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
    num_iters: u32,
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

#[cfg(feature = "validate")]
pub fn confirm_volume_exclusion(
    vs: &[V2D; NVERTS],
    contacts: &[ContactData],
    msg: &str,
) -> Result<(), String> {
    use crate::math::v2d::poly_to_string;
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

pub fn enforce_volume_exclusion(
    old_vs: &[V2D; NVERTS],
    mut new_vs: [V2D; NVERTS],
    contacts: Vec<ContactData>,
) -> Result<[V2D; NVERTS], String> {
    #[cfg(feature = "validate")]
    confirm_volume_exclusion(&old_vs, &contacts, "old_vs")?;

    for vi in 0..NVERTS {
        let old_v = old_vs[vi];
        let new_v = new_vs[vi];
        let new_u = new_vs[circ_ix_minus(vi, NVERTS)];
        let new_w = new_vs[circ_ix_plus(vi, NVERTS)];
        new_vs[vi] = move_point_out(
            &new_u, new_v, &new_w, old_v, &contacts, 20,
        );
    }

    #[cfg(feature = "validate")]
    confirm_volume_exclusion(&new_vs, &contacts, "new_vs")?;

    Ok(new_vs)
}

impl Cell {
    pub fn new(
        ix: usize,
        group_ix: usize,
        core: Core,
        interactions: &Interactions,
        parameters: &Parameters,
        rng: &mut Pcg32,
    ) -> Cell {
        let geom = core.calc_geom_state();
        let mech = core.calc_mech_state(&geom, parameters);
        let rac_rand = if parameters.randomization {
            RacRandState::new(rng, parameters)
        } else {
            RacRandState::default()
        };
        let chem = core.calc_chem_state(
            &mech,
            &rac_rand,
            &interactions,
            parameters,
        );
        Cell {
            ix,
            group_ix,
            core,
            rac_rand,
            geom,
            chem,
            mech,
        }
    }

    #[allow(unused)]
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
            state = state + dt * delta;
        }
        let geom_state = state.calc_geom_state();
        let mech_state =
            state.calc_mech_state(&geom_state, parameters);
        let chem_state = state.calc_chem_state(
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

        #[cfg(feature = "validate")]
        state.validate("euler", &parameters)?;

        Ok(Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(tstep, rng, parameters),
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
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
            max_iters: 100,
            atol: 1e-8,
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

        #[cfg(feature = "validate")]
        state.validate("rkdp5", &parameters)?;

        Ok(Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            core: state,
            rac_rand: self.rac_rand.update(tstep, rng, parameters),
            geom: geom_state,
            chem: chem_state,
            mech: mech_state,
        })
    }
}

/// Calculate the area of an "ideal" initial cell of radius R, if it has n vertices.
pub fn calc_init_cell_area(r: f64) -> f64 {
    let poly_coords = (0..NVERTS)
        .map(|vix| {
            let theta = (vix as f64) / (NVERTS as f64) * 2.0 * PI;
            V2D {
                x: r * theta.cos(),
                y: r * theta.sin(),
            }
        })
        .collect::<Vec<V2D>>();
    calc_poly_area(&poly_coords)
}
