// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::cell::chemistry::RacRandState;

use crate::cell::states::Core;
use crate::hardio::py_compare::{IntStepData, Writer};
use crate::interactions::{ContactData, Interactions};
use crate::math::geometry::calc_poly_area;

use crate::math::close_to_zero;
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::display::stringify_f64_arr;
use crate::utils::pcg32::Pcg32;

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
    /// State of Random Rac1 activity that affected `core`.
    pub rac_rand: RacRandState,
    /// Core state of the cell (position, Rho GTPase).
    pub core: Core,
    old_x_coas: [f64; NVERTS],
    old_x_cils: [f64; NVERTS],
    print_opts: PrintOptions,
}

#[derive(
    Copy, Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct PrintOptions {
    pub deltas: bool,
    pub interaction_updates: bool,
}

impl PrintOptions {
    pub fn any(&self) -> bool {
        self.deltas || self.interaction_updates
    }
}

impl Cell {
    pub fn new(
        ix: usize,
        group_ix: usize,
        core: Core,
        interactions: &Interactions,
        parameters: &Parameters,
        rng: &mut Pcg32,
        print_opts: PrintOptions,
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
            old_x_coas: interactions.x_coas,
            old_x_cils: interactions.x_cils,
            print_opts,
        }
    }

    pub fn print_tstep_header(&self, tstep: u32) {
        if self.print_opts.any() {
            println!("-----------------------------------");
            println!("tstep: {}, cell: {}", tstep, self.ix);
        }
    }

    pub fn print_interaction_updates(
        &self,
        updates: &[bool; NVERTS],
        old_interact_factors: &[f64; NVERTS],
        new_interact_factors: &[f64; NVERTS],
        description: &str,
    ) {
        if self.print_opts.interaction_updates {
            if updates.iter().any(|&x| x) {
                println!(
                    "old_coa = {}",
                    stringify_f64_arr(old_interact_factors, 4)
                );
                println!(
                    "new_coa = {}",
                    stringify_f64_arr(new_interact_factors, 4)
                );
            } else {
                println!("{}: no change", description);
                println!("{}: no change", description);
            }
        }
    }

    pub fn find_updates(
        old_interacts: &[f64; NVERTS],
        new_interacts: &[f64; NVERTS],
    ) -> [bool; NVERTS] {
        let mut updates = [false; NVERTS];
        for i in 0..NVERTS {
            updates[i] =
                !close_to_zero(old_interacts[i] - new_interacts[i]);
        }
        updates
    }

    #[allow(clippy::print_with_newline)]
    pub fn print_poly_delta_header(&self, old_state: &Core) {
        if self.print_opts.deltas {
            print!("-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
            print!(
                "init poly[0]: {:?}\n",
                [old_state.poly[0].x, old_state.poly[0].y]
            );
            print!(
                "init poly[15]: {:?}\n",
                [old_state.poly[15].x, old_state.poly[15].y]
            );
        }
    }

    pub fn write_state(
        &self,
        interactions: &Interactions,
        parameters: &Parameters,
        cil_updates: [bool; NVERTS],
        coa_updates: [bool; NVERTS],
        writer: &mut Writer,
    ) {
        let geom_state = self.core.geom;
        let mech_state = self.core.calc_mech_state(parameters);
        let chem_state = self.core.calc_chem_state(
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        let save_data = IntStepData {
            poly: self
                .core
                .poly
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            rac_acts: self.core.rac_acts,
            rac_inacts: self.core.rac_inacts,
            rho_acts: self.core.rho_acts,
            rho_inacts: self.core.rho_inacts,
            sum_forces: mech_state
                .sum_forces
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            uivs: geom_state
                .unit_in_vecs
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            kgtps_rac: chem_state.kgtps_rac,
            kdgtps_rac: chem_state.kdgtps_rac,
            kgtps_rho: chem_state.kgtps_rho,
            kdgtps_rho: chem_state.kdgtps_rho,
            rgtp_forces: mech_state
                .rgtp_forces
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            edge_forces: mech_state
                .edge_forces
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            uevs: geom_state
                .unit_edge_vecs
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            cyto_forces: mech_state
                .cyto_forces
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            x_cils: interactions.x_cils,
            x_coas: interactions.x_coas,
            rac_act_net_fluxes: chem_state.rac_act_net_fluxes,
            edge_strains: mech_state.edge_strains,
            poly_area: calc_poly_area(&self.core.poly),
            coa_updates,
            cil_updates,
        };
        writer.save_int_step(save_data);
    }

    /// Executes the same logic as `Cell::simulate_euler`, but saves
    /// additional data for comparison with the Python model.
    #[allow(clippy::too_many_arguments)]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        num_int_steps: u32,
        interactions: &Interactions,
        contacts: Vec<ContactData>,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
        rng: &mut Pcg32,
        writer: &mut Writer,
    ) -> Result<Cell, String> {
        self.print_tstep_header(tstep);

        let coa_updates = Self::find_updates(
            &self.old_x_coas,
            &interactions.x_coas,
        );
        self.print_interaction_updates(
            &coa_updates,
            &self.old_x_coas,
            &interactions.x_coas,
            "coa",
        );

        let cil_updates = Self::find_updates(
            &self.old_x_cils,
            &interactions.x_cils,
        );

        self.print_interaction_updates(
            &cil_updates,
            &self.old_x_cils,
            &interactions.x_cils,
            "cil",
        );

        let dt = 1.0 / (num_int_steps as f64);
        let mut state = self.core;
        for _ in 0..num_int_steps {
            self.print_poly_delta_header(&state);
            self.write_state(
                interactions,
                parameters,
                cil_updates,
                coa_updates,
                writer,
            );
            // d(state)/dt = dynamics_f(state) <- calculate RHS of ODE
            let delta = state.derivative(
                &self.rac_rand,
                interactions,
                world_parameters,
                parameters,
            );
            let mut new_state = state + delta.time_step(dt);
            let delta_poly_0 = new_state.poly[0] - state.poly[0];
            let delta_poly_15 = new_state.poly[15] - state.poly[15];
            // Enforcing volume exclusion! Tricky!
            new_state
                .enforce_volume_exclusion(&state.poly, &contacts)?;
            if self.print_opts.deltas {
                println!(
                    "calculated Delta poly(0): {}",
                    delta_poly_0
                );
                println!(
                    "Delta poly after VE: {}",
                    new_state.poly[0] - state.poly[0]
                );
                println!("final poly[0]: {:?}", new_state.poly[0]);
                println!("actual Delta poly(15): {}", delta_poly_15);
                println!(
                    "Delta poly after VE: {}",
                    new_state.poly[15] - state.poly[15]
                );
                println!("final poly[15]: {:?}", new_state.poly[15]);
            }
            state = new_state;
        }

        #[cfg(feature = "validate")]
        state.validate("euler")?;

        Ok(Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            rac_rand: self.rac_rand.update(
                tstep + 1,
                rng,
                parameters,
            ),
            core: state,
            old_x_coas: interactions.x_coas,
            old_x_cils: interactions.x_cils,
            print_opts: self.print_opts,
        })
    }
}
