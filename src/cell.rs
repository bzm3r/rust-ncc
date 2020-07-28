// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::chemistry::{
    calc_conc_rgtps, calc_kdgtps_rac, calc_kdgtps_rho,
    calc_kgtps_rac, calc_kgtps_rho, calc_net_fluxes, gen_rgtp_distrib, RgtpLayout,
};
use crate::consts::NVERTS;
use crate::math::{hill_function3, max_f32, min_f32, P2D};
use crate::mechanics::{calc_cyto_forces, calc_edge_forces, calc_edge_vecs, calc_rgtp_forces};
use crate::parameters::Parameters;
use crate::rkdp5::{rkdp5, AuxArgs};
use crate::utils::circ_ix_minus;
use crate::world::InteractionState;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;
use std::ops::{Add, Div, Mul, Sub};
use std::fmt::Display;
use std::fmt;

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct CellState {
    pub vertex_coords: [P2D; NVERTS as usize],
    rac_acts: [f32; NVERTS as usize],
    rac_inacts: [f32; NVERTS as usize],
    rho_acts: [f32; NVERTS as usize],
    rho_inacts: [f32; NVERTS as usize],
}

impl Add for CellState {
    type Output = CellState;

    fn add(self, rhs: CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = &self.vertex_coords[i] + &rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] + &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] + &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] + &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] + &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Sub for CellState {
    type Output = CellState;

    fn sub(self, rhs: CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = &self.vertex_coords[i] - &rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] - &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] - &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] - &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] - &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Add for &CellState {
    type Output = CellState;

    fn add(self, rhs: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = &self.vertex_coords[i] + &rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] + &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] + &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] + &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] + &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Add<&CellState> for CellState {
    type Output = CellState;

    fn add(self, rhs: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = &self.vertex_coords[i] + &rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] + &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] + &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] + &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] + &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Add<CellState> for &CellState {
    type Output = CellState;

    fn add(self, rhs: CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = &self.vertex_coords[i] + &rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] + &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] + &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] + &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] + &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Div for CellState {
    type Output = CellState;

    fn div(self, rhs: CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = self.vertex_coords[i] / rhs.vertex_coords[i];
            rac_acts[i] = &self.rac_acts[i] / &rhs.rac_acts[i];
            rac_inacts[i] = &self.rac_inacts[i] / &rhs.rac_inacts[i];
            rho_acts[i] = &self.rho_acts[i] / &rhs.rho_acts[i];
            rho_inacts[i] = &self.rho_inacts[i] / &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Mul<CellState> for f32 {
    type Output = CellState;

    fn mul(self, rhs: CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = self * rhs.vertex_coords[i];
            rac_acts[i] = self * &rhs.rac_acts[i];
            rac_inacts[i] = self * &rhs.rac_inacts[i];
            rho_acts[i] = self * &rhs.rho_acts[i];
            rho_inacts[i] = self * &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Mul<&CellState> for f32 {
    type Output = CellState;

    fn mul(self, rhs: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = self * rhs.vertex_coords[i];
            rac_acts[i] = self * &rhs.rac_acts[i];
            rac_inacts[i] = self * &rhs.rac_inacts[i];
            rho_acts[i] = self * &rhs.rho_acts[i];
            rho_inacts[i] = self * &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

impl Mul<&CellState> for &f32 {
    type Output = CellState;

    fn mul(self, rhs: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = self * rhs.vertex_coords[i];
            rac_acts[i] = self * &rhs.rac_acts[i];
            rac_inacts[i] = self * &rhs.rac_inacts[i];
            rho_acts[i] = self * &rhs.rho_acts[i];
            rho_inacts[i] = self * &rhs.rho_inacts[i]
        }

        Self::Output {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct RacRandomState {
    x_rands: [f32; NVERTS as usize],
    next_updates: [u32; NVERTS as usize],
}

impl RacRandomState {
    fn init() -> RacRandomState {
        RacRandomState::default()
    }

    fn update(&self, _tstep: u32) -> RacRandomState {
        RacRandomState::default()
    }
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct MechState {
    pub edge_strains: [f32; NVERTS as usize],
    pub rgtp_forces: [P2D; NVERTS as usize],
    pub cyto_forces: [P2D; NVERTS as usize],
    pub edge_forces: [P2D; NVERTS as usize],
    pub avg_tens_strain: f32
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct ChemState {
    pub kdgtps_rac: [f32; NVERTS as usize],
    pub kgtps_rac: [f32; NVERTS as usize],
    pub rac_act_net_fluxes: [f32; NVERTS as usize],
    pub rac_inact_net_fluxes: [f32; NVERTS as usize],
    pub kdgtps_rho: [f32; NVERTS as usize],
    pub kgtps_rho: [f32; NVERTS as usize],
    pub rac_cyto: f32,
    pub rho_cyto: f32,
    pub rho_act_net_fluxes: [f32; NVERTS as usize],
    pub rho_inact_net_fluxes: [f32; NVERTS as usize],
    pub x_tens: f32
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct GeomState {
    pub unit_edge_vecs: [P2D; NVERTS as usize],
    pub edge_lens: [f32; NVERTS as usize],
    pub unit_inward_vecs: [P2D; NVERTS as usize],
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct CellDepVars {
    geom_state: GeomState,
    chem_state: ChemState,
    mech_state: MechState,
}

fn fmt_var_arr<T: fmt::Display>(f: &mut fmt::Formatter<'_>, description: &str, vars: &[T; NVERTS as usize]) -> fmt::Result {
    let contents = vars.iter().map(|x| format!("{}", x)).collect::<Vec<String>>().join(", ");
    write!(f, "{}: [{}]\n", description, contents)
}

impl Display for CellState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        println!("----------");
        fmt_var_arr(f, "vertex_coords", &self.vertex_coords)?;
        fmt_var_arr(f, "rac_acts", &self.rac_acts)?;
        fmt_var_arr(f, "rac_inacts", &self.rac_inacts)?;
        fmt_var_arr(f, "rho_acts", &self.rho_acts)?;
        fmt_var_arr(f, "rho_inacts", &self.rho_inacts)

    }
}

impl Display for GeomState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "edge_lens", &self.edge_lens)?;
        fmt_var_arr(f, "uivs", &self.unit_inward_vecs)?;
        fmt_var_arr(f, "uevs", &self.unit_edge_vecs)
    }
}

impl Display for MechState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "rgtp_forces", &self.rgtp_forces)?;
        fmt_var_arr(f, "edge_strains", &self.edge_strains)?;
        write!(f, "avg_tens_strain: {}\n", self.avg_tens_strain)?;
        fmt_var_arr(f, "edge_forces", &self.edge_forces)?;
        fmt_var_arr(f, "cyto_forces", &self.cyto_forces)
    }
}


impl Display for ChemState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "x_tens: {}\n", self.x_tens)?;
        fmt_var_arr(f, "kgtps_rac", &self.kgtps_rac)?;
        fmt_var_arr(f, "kdgtps_rac", &self.kdgtps_rac)?;
        fmt_var_arr(f, "kgtps_rho", &self.kgtps_rho)?;
        fmt_var_arr(f, "kdgtps_rho", &self.kdgtps_rac)
    }
}

impl Display for CellDepVars {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\n", &self.geom_state)?;
        write!(f, "{}\n", &self.mech_state)?;
        write!(f, "{}\n", &self.chem_state)
    }
}

impl CellState {
    pub fn num_vars() -> usize {
        (NVERTS * 6) as usize
    }

    pub(crate) fn calc_geom_state(state: &CellState) -> GeomState {
        let evs = calc_edge_vecs(&state.vertex_coords);
        let mut edge_lens = [0.0_f32; NVERTS as usize];
        (0..NVERTS as usize).for_each(|i| edge_lens[i] = (&evs[i]).mag());
        let mut uevs = [P2D::default(); NVERTS as usize];
        (0..NVERTS as usize).for_each(|i| uevs[i] = (&evs[i]).unitize());
        let mut uivs = [P2D::default(); NVERTS as usize];
        (0..NVERTS as usize).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS as usize);
            let tangent = (uevs[i] + uevs[im1]).unitize();
            uivs[i] = tangent.normal();
        });

        GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        }
    }

    pub fn calc_mech_state(
        state: &CellState,
        geom_state: &GeomState,
        parameters: &Parameters,
    ) -> MechState {
        let GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        } = geom_state;
        let rgtp_forces = calc_rgtp_forces(
            &state.rac_acts,
            &state.rho_acts,
            uivs,
            parameters.halfmax_vertex_rgtp_act,
            parameters.const_protrusive,
            parameters.const_retractive,
        );
        let cyto_forces = calc_cyto_forces(
            &state.vertex_coords,
            &uivs,
            parameters.rest_area,
            parameters.stiffness_ctyo,
        );
        let mut edge_strains = [0.0_f32; NVERTS as usize];
        (0..NVERTS as usize)
            .for_each(|i| edge_strains[i] = (edge_lens[i] / parameters.rest_edge_len) - 1.0);
        let edge_forces = calc_edge_forces(&edge_strains, uevs, parameters.stiffness_edge);
        let avg_tens_strain = edge_strains.iter().map(|&es| if es < 0.0 {
            0.0
        } else {
            es
        }).sum::<f32>() / NVERTS as f32;
        MechState {
            edge_strains,
            rgtp_forces,
            cyto_forces,
            edge_forces,
            avg_tens_strain,
        }
    }

    pub fn calc_chem_state(
        state: &CellState,
        geom_state: &GeomState,
        mech_state: &MechState,
        rac_rand_state: &RacRandomState,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> ChemState {
        let GeomState { edge_lens, .. } = geom_state;
        let mut avg_edge_lens: [f32; NVERTS as usize] = [0.0_f32; NVERTS as usize];
        (0..NVERTS as usize).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS as usize);
            avg_edge_lens[i] = (edge_lens[i] + edge_lens[im1]) / 2.0;
        });

        let conc_rac_acts = calc_conc_rgtps(&avg_edge_lens, &state.rac_acts);
        let conc_rac_inacts = calc_conc_rgtps(&avg_edge_lens, &state.rac_inacts);
        let conc_rho_acts = calc_conc_rgtps(&avg_edge_lens, &state.rho_acts);
        let conc_rho_inacts = calc_conc_rgtps(&avg_edge_lens, &state.rho_inacts);

        let kgtps_rac = calc_kgtps_rac(
            &state.rac_acts,
            &conc_rac_acts,
            &rac_rand_state.x_rands,
            &inter_state.x_coas,
            &inter_state.x_chemoas,
            parameters.kgtp_rac,
            parameters.kgtp_rac_auto,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let MechState { avg_tens_strain, .. } = mech_state;
        let x_tens = parameters.tension_inhib * hill_function3(
            parameters.halfmax_tension_inhib,
            *avg_tens_strain,
        );
        let kdgtps_rac = calc_kdgtps_rac(
            &state.rac_acts,
            &conc_rho_acts,
            &inter_state.x_cils,
            x_tens,
            parameters.kdgtp_rac,
            parameters.kdgtp_rho_on_rac,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let kgtps_rho = calc_kgtps_rho(
            &state.rho_acts,
            &conc_rho_acts,
            &inter_state.x_cils,
            parameters.kgtp_rho,
            parameters.halfmax_vertex_rgtp_conc,
            parameters.kgtp_rho_auto,
        );
        let kdgtps_rho = calc_kdgtps_rho(
            &state.rho_acts,
            &conc_rac_acts,
            parameters.kdgtp_rho,
            parameters.kdgtp_rac_on_rho,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let rac_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rac_acts,
        );
        let rho_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rho_acts,
        );
        let rac_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rac_inacts,
        );
        let rho_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            parameters.diffusion_rgtp,
            &conc_rho_inacts,
        );

        let rac_cyto =
            1.0 - state.rac_acts.iter().sum::<f32>() - state.rac_inacts.iter().sum::<f32>();
        let rho_cyto =
            1.0 - state.rho_acts.iter().sum::<f32>() - state.rho_inacts.iter().sum::<f32>();
        ChemState {
            x_tens,
            kdgtps_rac,
            kgtps_rac,
            rac_act_net_fluxes,
            rac_inact_net_fluxes,
            rho_act_net_fluxes,
            rho_inact_net_fluxes,
            kdgtps_rho,
            kgtps_rho,
            rac_cyto,
            rho_cyto,
        }
    }

    fn calc_dep_vars(
        state: &CellState,
        rac_rand_state: &RacRandomState,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> CellDepVars {
        let geom_state = Self::calc_geom_state(state);
        let mech_state = Self::calc_mech_state(state, &geom_state, parameters);
        let chem_state = Self::calc_chem_state(
            state,
            &geom_state,
            &mech_state,
            rac_rand_state,
            inter_state,
            parameters,
        );

        CellDepVars {
            geom_state,
            chem_state,
            mech_state,
        }
    }

    pub fn dynamics_f(
        state: &CellState,
        rac_rand_state: &RacRandomState,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> CellState {
        let CellDepVars {
            chem_state,
            mech_state,
            ..
        } = CellState::calc_dep_vars(state, rac_rand_state, inter_state, parameters);
        let mut delta = CellState::default();
        for i in 0..NVERTS as usize {
            let inactivated_rac = chem_state.kdgtps_rac[i] * state.rac_acts[i];
            let activated_rac = chem_state.kgtps_rac[i] * state.rac_inacts[i];
            let delta_rac_activated = activated_rac - inactivated_rac;
            let rac_cyto_exchange = {
                let rac_mem_on = parameters.k_mem_on_vertex * chem_state.rac_cyto;
                let rac_mem_off = parameters.k_mem_off * state.rac_inacts[i];
                rac_mem_on - rac_mem_off
            };
            let vertex_rac_act_flux = chem_state.rac_act_net_fluxes[i];
            let vertex_rac_inact_flux = chem_state.rac_inact_net_fluxes[i];
            delta.rac_acts[i] = delta_rac_activated + vertex_rac_act_flux;
            delta.rac_inacts[i] = rac_cyto_exchange + vertex_rac_inact_flux - delta_rac_activated;

            let inactivated_rho = chem_state.kdgtps_rho[i] * state.rho_acts[i];
            let activated_rho = chem_state.kgtps_rho[i] * state.rho_inacts[i];
            let delta_rho_activated = activated_rho - inactivated_rho;
            let rho_cyto_exchange = {
                let rho_mem_on = parameters.k_mem_on_vertex * chem_state.rho_cyto;
                let rho_mem_off = parameters.k_mem_off * state.rho_inacts[i];
                rho_mem_on - rho_mem_off
            };
            let vertex_rho_act_flux = chem_state.rho_act_net_fluxes[i];
            let vertex_rho_inact_flux = chem_state.rho_inact_net_fluxes[i];
            delta.rho_acts[i] = delta_rho_activated + vertex_rho_act_flux;
            delta.rho_inacts[i] = rho_cyto_exchange + vertex_rho_inact_flux - delta_rho_activated;

            let sum_f =
                mech_state.rgtp_forces[i] + mech_state.cyto_forces[i] + mech_state.edge_forces[i]
                    - mech_state.edge_forces[circ_ix_minus(i as usize, NVERTS as usize)];
            let delta_x = (1.0 / parameters.vertex_eta) * sum_f.x;
            delta.vertex_coords[i] = (1.0 / parameters.vertex_eta) * sum_f;
        }
        delta
    }

    // pub fn decompress_force_data(
    //     &self,
    //     rgtp_forces: bool,
    //     edge_forces: bool,
    //     cyto_forces: bool,
    //     parameters: &Parameters,
    // ) -> ForceData {
    //     let (rgtp_fs, edge_fs, cyto_fs) = if rgtp_forces || edge_forces {
    //         let pg = self.calc_geom_vars();
    //         let rfs = if rgtp_forces {
    //             Some(calc_rgtp_forces(
    //                 &self.rac_acts,
    //                 &self.rho_acts,
    //                 &pg.unit_inward_vecs,
    //                 parameters.halfmax_vertex_rgtp_act,
    //                 parameters.const_protrusive,
    //                 parameters.const_retractive,
    //             ))
    //         } else {
    //             None
    //         };
    //
    //         let efs = if edge_forces {
    //             let edge_strains =
    //                 calc_edge_strains(&pg.edge_lens, parameters.rest_edge_len);
    //             Some(calc_edge_forces(
    //                 &edge_strains,
    //                 &pg.unit_edge_vecs,
    //                 parameters.stiffness_edge,
    //             ))
    //         } else {
    //             None
    //         };
    //
    //         let cfs = if cyto_forces {
    //             Some(calc_cyto_forces(
    //                 &self.vertex_coords,
    //                 &pg.unit_inward_vecs,
    //                 parameters.rest_area,
    //                 parameters.stiffness_ctyo,
    //             ))
    //         } else {
    //             None
    //         };
    //
    //         (rfs, efs, cfs)
    //     } else {
    //         (None, None, None)
    //     };
    //
    //     ForceData {
    //         rgtp_fs,
    //         edge_fs,
    //         cyto_fs,
    //     }
    // }

    pub fn new(
        vertex_coords: [P2D; NVERTS as usize],
        rac_acts: [f32; NVERTS as usize],
        rac_inacts: [f32; NVERTS as usize],
        rho_acts: [f32; NVERTS as usize],
        rho_inacts: [f32; NVERTS as usize],
    ) -> CellState {
        // x_cils: [f32; NVERTS], x_coas: [f32; NVERTS], x_chemoas: [f32; NVERTS], x_rands: [f32; NVERTS], x_bdrys: [f32; NVERTS as usize];
        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn scalar_mul(&self, s: f32) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = s * self.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] * s;
            rac_inacts[i] = self.rac_inacts[i] * s;
            rho_acts[i] = self.rho_acts[i] * s;
            rho_inacts[i] = self.rho_inacts[i] * s;
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn scalar_add(&self, s: f32) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = s + self.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] + s;
            rac_inacts[i] = self.rac_inacts[i] + s;
            rho_acts[i] = self.rho_acts[i] + s;
            rho_inacts[i] = self.rho_inacts[i] + s;
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn abs(&self) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = vertex_coords[i].abs();
            rac_acts[i] = self.rac_acts[i].abs();
            rac_inacts[i] = self.rac_inacts[i].abs();
            rho_acts[i] = self.rho_acts[i].abs();
            rho_inacts[i] = self.rho_inacts[i].abs();
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn powi(&self, x: i32) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = vertex_coords[i].powi(x);
            rac_acts[i] = self.rac_acts[i].powi(x);
            rac_inacts[i] = self.rac_inacts[i].powi(x);
            rho_acts[i] = self.rho_acts[i].powi(x);
            rho_inacts[i] = self.rho_inacts[i].powi(x);
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn max(&self, other: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = vertex_coords[i].max(&other.vertex_coords[i]);
            rac_acts[i] = max_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] = max_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] = max_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] = max_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn min(&self, other: &CellState) -> CellState {
        let mut vertex_coords = [P2D::default(); NVERTS as usize];
        let mut rac_acts = [0.0_f32; NVERTS as usize];
        let mut rac_inacts = [0.0_f32; NVERTS as usize];
        let mut rho_acts = [0.0_f32; NVERTS as usize];
        let mut rho_inacts = [0.0_f32; NVERTS as usize];

        for i in 0..(NVERTS as usize) {
            vertex_coords[i] = vertex_coords[i].min(&other.vertex_coords[i]);
            rac_acts[i] = min_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] = min_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] = min_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] = min_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn sum(&self) -> f32 {
        let mut r: f32 = 0.0;

        for i in 0..(NVERTS as usize) {
            r += self.vertex_coords[i].x + self.vertex_coords[i].y;
            r += self.rac_acts[i];
            r += self.rac_inacts[i];
            r += self.rho_acts[i];
            r += self.rho_inacts[i];
        }

        r
    }

    pub fn average(&self) -> f32 {
        self.sum() / (Self::num_vars() as f32)
    }

    pub fn validate(&self, loc_str: &str) {
        if self.rac_acts.iter().any(|&r| r < 0.0_f32) {
            panic!("{}: neg rac_acts: {:?}", loc_str, self.rac_acts);
        }
        if self.rac_inacts.iter().any(|&r| r < 0.0_f32) {
            panic!("{}: neg rac_inacts: {:?}", loc_str, self.rac_inacts);
        }
        if self.rho_acts.iter().any(|&r| r < 0.0_f32) {
            panic!("{}: neg rho_acts: {:?}", loc_str, self.rho_acts);
        }
        if self.rho_inacts.iter().any(|&r| r < 0.0_f32) {
            panic!("{}: neg rho_inacts: {:?}", loc_str, self.rho_inacts);
        }
        let sum_rac_mem = self.rac_inacts.iter().sum::<f32>() + self.rac_acts.iter().sum::<f32>();
        if sum_rac_mem > 1.0 || sum_rac_mem < 0.0 {
            panic!("{}: problem in sum of rac_mem: {}", loc_str, sum_rac_mem);
        }
        let sum_rho_mem = self.rho_inacts.iter().sum::<f32>() + self.rho_acts.iter().sum::<f32>();
        if sum_rho_mem > 1.0 || sum_rho_mem < 0.0 {
            panic!("{}: problem in sum of rho_mem: {}", loc_str, sum_rho_mem);
        }
        println!("{}: successfully validated", loc_str)
    }
}

/// Model cell.
#[derive(Clone, Deserialize, Serialize, Schematize)]
pub struct Cell {
    /// Cell index.
    pub ix: u32,
    /// Index of cell type.
    pub group_ix: u32,
    // /// If `true`, `evolve` does not change current cell state.
    // skip_dynamics: bool,
    pub state: CellState,
    pub rac_randomization: RacRandomState,
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
    ) -> Cell {
        let vertex_coords = match vertex_gen_info {
            VertexGenInfo::Centroid(centroid) => gen_vertex_coords(centroid, parameters.cell_r),
            VertexGenInfo::Coordinates(vcs) => vcs,
        };
        let (rac_acts, rac_inacts) = gen_rgtp_distrib(
            parameters.init_frac_active,
            parameters.init_frac_inactive,
            //&RgtpLayout::Random,
            &RgtpLayout::BiasedVertices(vec![0, 1, 2, 3]),
        );
        let (rho_acts, rho_inacts) = gen_rgtp_distrib(
            parameters.init_frac_active,
            parameters.init_frac_inactive,
            //&RgtpLayout::Random,
            &RgtpLayout::BiasedVertices(vec![8, 9, 10, 11]),
        );
        let state = CellState::new(vertex_coords, rac_acts, rac_inacts, rho_acts, rho_inacts);
        let rac_rand_state = RacRandomState::init();
        Cell {
            ix,
            group_ix,
            state,
            rac_randomization: rac_rand_state,
        }
    }

    #[allow(unused)]
    pub fn simulate_euler(
        &self,
        tstep: u32,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> Cell {
        let mut state = self.state.clone();
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for i in 0..nsteps {
            println!("++++++++++++");
            println!("{}", state);
            let dep_vars = CellState::calc_dep_vars(&state, &self.rac_randomization, inter_state, parameters);
            println!("{}", dep_vars);
            println!("++++++++++++");
            let delta = CellState::dynamics_f(
                &state,
                &self.rac_randomization,
                &inter_state,
                parameters,
            );
            state = state + dt * delta;
        }
        // println!("++++++++++++");
        // #[cfg(debug_assertions)]
        // state.validate("euler");
        // println!("{}", state);
        // let dep_vars = CellState::calc_dep_vars(&state, &self.rac_randomization, inter_state, parameters);
        // println!("{}", dep_vars);
        // println!("++++++++++++");
        Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_randomization: self.rac_randomization.update(tstep),
        }
    }

    pub fn simulate_rkdp5(
        &self,
        tstep: u32,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> Cell {
        println!("using rkdp5...");
        let aux_args = AuxArgs {
            max_iters: 100,
            atol: 1e-8,
            rtol: 1e-3,
            init_h_factor: Some(0.1),
        };
        let result = rkdp5(
            1.0,
            CellState::dynamics_f,
            &self.state,
            &self.rac_randomization,
            inter_state,
            parameters,
            aux_args,
        );

        println!(
            "num_iters: {}, num_rejections: {}",
            result.num_iters, result.num_rejections
        );
        let state = result.y.expect("too many iterations!");
        #[cfg(debug_assertions)]
        state.validate("rkdp5");
        println!("{}", state);
        let dep_vars = CellState::calc_dep_vars(&state, &self.rac_randomization, inter_state, parameters);
        println!("{}", dep_vars);
        Cell {
            ix: self.ix,
            group_ix: self.group_ix,
            state,
            rac_randomization: self.rac_randomization.update(tstep),
        }
    }
}
