// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::cell::chemistry::RacRandState;

use crate::cell::states::{Core, GeomState, MechState};
use crate::cell::Cell;
use crate::hardio::pycomp::{IntStepData, Writer};
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
pub struct PyCompCell {
    cell: Cell,
    old_x_coas: [f64; NVERTS],
    old_x_cils: [f64; NVERTS],
    print_opts: PrintOptions,
}

#[derive(
    Copy, Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct PrintOptions {
    deltas: bool,
    interaction_updates: bool,
}

impl PrintOptions {
    pub fn any(&self) -> bool {
        self.deltas || self.interaction_updates
    }
}

impl PyCompCell {
    pub fn new(
        ix: usize,
        group_ix: usize,
        core: Core,
        interactions: &Interactions,
        parameters: &Parameters,
        rng: &mut Pcg32,
        print_opts: PrintOptions,
    ) -> PyCompCell {
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
        PyCompCell {
            cell: Cell {
                ix,
                group_ix,
                core,
                rac_rand,
                geom,
                chem,
                mech,
            },
            old_x_coas: interactions.x_coas,
            old_x_cils: interactions.x_cils,
            print_opts,
        }
    }

    pub fn print_tstep_header(&self, tstep: u32) {
        if self.print_opts.any() {
            println!("-----------------------------------");
            println!("tstep: {}, cell: {}", tstep, self.cell.ix);
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
            updates[i] = !close_to_zero(
                old_interacts[i] - new_interacts.[i],
            );
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

    pub fn write_state(&self, parameters: &Parameters, writer: &mut Writer) {
        let geom_state = self.cell.core.calc_geom_state();
        let mech_state = self.cell.core.calc_mech_state(parameters);
        let chem_state = self.cell.core.calc_chem_state(
            &mech_state,
            &self.rac_rand,
            &interactions,
            parameters,
        );
        let save_data = IntStepData {
            poly: self
                .cell
                .state
                .poly
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            rac_acts: self.cell.state.rac_acts,
            rac_inacts: self.cell.state.rac_inacts,
            rho_acts: self.cell.state.rho_acts,
            rho_inacts: self.cell.state.rho_inacts,
            sum_forces: mech_state
                .sum_forces
                .iter()
                .map(|v| [v.x, v.y])
                .collect::<Vec<[f64; 2]>>(),
            uivs: geom_state
                .unit_inward_vecs
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
            edge_forces_minus: mech_state
                .edge_forces_minus
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
            conc_rac_acts: chem_state.conc_rac_acts,
            x_cils: interactions.x_cils,
            x_coas: interactions.x_coas,
            rac_act_net_fluxes: chem_state.rac_act_net_fluxes,
            edge_strains: mech_state.edge_strains,
            poly_area: calc_poly_area(&state.poly),
            coa_updates,
            cil_update,
        };
        writer.save_int_step(save_data);
    }

    #[allow(unused)]
    /// A wrapper around `Cell::simulate_euler`, that outputs data for
    /// comparison with the Python model.
    pub fn simulate_euler(
        &self,
        tstep: u32,
        num_int_steps: u32,
        interactions: &Interactions,
        contact_data: Vec<ContactData>,
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

        let mut cil_updates = Self::find_updates(
            &self.old_x_cils,
            &interactions.x_cils,
        );
        self.print_interaction_updates(
            &cil_updates,
            &self.old_x_cils,
            &interactions.x_cils,
        );

        for int_step in 0..num_int_steps {
            self.print_poly_delta_header(&self.cell.state);
            self.write_state(writer);
            // d(state)/dt = dynamics_f(state) <- calculate RHS of ODE
            let delta = Core::derivative();
            state = state + dt * delta;
            let actual_dp_0 = state.poly[0] - self.core.poly[0];
            let actual_dp_15 = state.poly[15] - self.core.poly[15];
            // Enforcing volume exclusion! Tricky!
            state.poly = enforce_volume_exclusion_py(
                talkative,
                &self.core.poly,
                state.poly,
                contact_data.clone(),
                dt,
                parameters.const_protrusive,
                world_parameters.vertex_eta,
                &geom_state.unit_inward_vecs,
            )?;
            if talkative {
                println!("actual Delta poly(0): {}", actual_dp_0);
                println!(
                    "Delta poly after VE: {}",
                    state.poly[0] - self.core.poly[0]
                );
                println!("final poly[0]: {:?}", state.poly[0]);
                println!("actual Delta poly(15): {}", actual_dp_15);
                println!(
                    "Delta poly after VE: {}",
                    state.poly[15] - self.core.poly[15]
                );
                println!("final poly[15]: {:?}", state.poly[15]);
            }
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
            old_x_coas: interactions.x_coas,
            old_x_cils: interactions.x_cils,
        })
    }
}
