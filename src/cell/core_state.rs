use crate::cell::chemistry::{
    calc_conc_rgtps, calc_kdgtps_rac, calc_kdgtps_rho, calc_kgtps_rac, calc_kgtps_rho,
    calc_net_fluxes, RacRandState, RgtpDistribution,
};
use crate::cell::mechanics::{
    calc_cyto_forces, calc_edge_forces, calc_edge_vecs, calc_rgtp_forces,
};
use crate::interactions::InteractionState;
use crate::math::p2d::P2D;
use crate::math::{hill_function3, max_f32, min_f32};
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::circ_ix_minus;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct CoreState {
    pub vertex_coords: [P2D; NVERTS],
    rac_acts: [f32; NVERTS],
    rac_inacts: [f32; NVERTS],
    rho_acts: [f32; NVERTS],
    rho_inacts: [f32; NVERTS],
}

impl Add for CoreState {
    type Output = CoreState;

    fn add(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.vertex_coords[i] + rhs.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] + rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] + rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] + rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] + rhs.rho_inacts[i]
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

impl Sub for CoreState {
    type Output = CoreState;

    fn sub(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.vertex_coords[i] - rhs.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] - rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] - rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] - rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] - rhs.rho_inacts[i]
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

impl Div for CoreState {
    type Output = CoreState;

    fn div(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self.vertex_coords[i] / rhs.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] / rhs.rac_acts[i];
            rac_inacts[i] = self.rac_inacts[i] / rhs.rac_inacts[i];
            rho_acts[i] = self.rho_acts[i] / rhs.rho_acts[i];
            rho_inacts[i] = self.rho_inacts[i] / rhs.rho_inacts[i]
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

impl Mul<CoreState> for f32 {
    type Output = CoreState;

    fn mul(self, rhs: CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = self * rhs.vertex_coords[i];
            rac_acts[i] = self * rhs.rac_acts[i];
            rac_inacts[i] = self * rhs.rac_inacts[i];
            rho_acts[i] = self * rhs.rho_acts[i];
            rho_inacts[i] = self * rhs.rho_inacts[i]
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
pub struct MechState {
    pub edge_strains: [f32; NVERTS],
    pub rgtp_forces: [P2D; NVERTS],
    pub cyto_forces: [P2D; NVERTS],
    pub edge_forces: [P2D; NVERTS],
    pub avg_tens_strain: f32,
    pub sum_fs: [P2D; NVERTS],
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct ChemState {
    pub kdgtps_rac: [f32; NVERTS],
    pub kgtps_rac: [f32; NVERTS],
    pub rac_act_net_fluxes: [f32; NVERTS],
    pub rac_inact_net_fluxes: [f32; NVERTS],
    pub kdgtps_rho: [f32; NVERTS],
    pub kgtps_rho: [f32; NVERTS],
    pub rac_cyto: f32,
    pub rho_cyto: f32,
    pub rho_act_net_fluxes: [f32; NVERTS],
    pub rho_inact_net_fluxes: [f32; NVERTS],
    pub x_tens: f32,
}

#[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
pub struct GeomState {
    pub unit_edge_vecs: [P2D; NVERTS],
    pub edge_lens: [f32; NVERTS],
    pub unit_inward_vecs: [P2D; NVERTS],
}

// #[derive(Copy, Clone, Debug, Default, Deserialize, Serialize, Schematize)]
// pub struct DepStates {
//     pub geom_state: GeomState,
//     pub chem_state: ChemState,
//     pub mech_state: MechState,
// }

pub fn fmt_var_arr<T: fmt::Display>(
    f: &mut fmt::Formatter<'_>,
    description: &str,
    vars: &[T; NVERTS],
) -> fmt::Result {
    let contents = vars
        .iter()
        .map(|x| format!("{}", x))
        .collect::<Vec<String>>()
        .join(", ");
    writeln!(f, "{}: [{}]", description, contents)
}

impl Display for CoreState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // println!("----------");
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
        writeln!(f, "avg_tens_strain: {}", self.avg_tens_strain)?;
        fmt_var_arr(f, "edge_forces", &self.edge_forces)?;
        fmt_var_arr(f, "cyto_forces", &self.cyto_forces)?;
        fmt_var_arr(f, "tot_forces", &self.sum_fs)
    }
}

impl Display for ChemState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "x_tens: {}", self.x_tens)?;
        fmt_var_arr(f, "kgtps_rac", &self.kgtps_rac)?;
        fmt_var_arr(f, "kdgtps_rac", &self.kdgtps_rac)?;
        fmt_var_arr(f, "kgtps_rho", &self.kgtps_rho)?;
        fmt_var_arr(f, "kdgtps_rho", &self.kdgtps_rac)
    }
}

// impl Display for DepStates {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         writeln!(f, "{}", &self.geom_state)?;
//         writeln!(f, "{}", &self.mech_state)?;
//         writeln!(f, "{}", &self.chem_state)
//     }
// }

impl CoreState {
    pub fn num_vars() -> usize {
        (NVERTS * 6) as usize
    }

    pub fn calc_geom_state(&self) -> GeomState {
        let evs = calc_edge_vecs(&self.vertex_coords);
        let mut edge_lens = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| edge_lens[i] = (&evs[i]).mag());
        let mut uevs = [P2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| uevs[i] = (&evs[i]).unitize());
        let mut uivs = [P2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS);
            let tangent = (uevs[i] + uevs[im1]).unitize();
            uivs[i] = tangent.normal();
        });

        GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        }
    }

    pub fn calc_mech_state(&self, geom_state: &GeomState, parameters: &Parameters) -> MechState {
        let GeomState {
            unit_edge_vecs: uevs,
            edge_lens,
            unit_inward_vecs: uivs,
        } = geom_state;
        let rgtp_forces = calc_rgtp_forces(
            &self.rac_acts,
            &self.rho_acts,
            uivs,
            parameters.halfmax_vertex_rgtp_act,
            parameters.const_protrusive,
            parameters.const_retractive,
        );
        let cyto_forces = calc_cyto_forces(
            &self.vertex_coords,
            &uivs,
            parameters.rest_area,
            parameters.stiffness_ctyo,
        );
        let mut edge_strains = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| edge_strains[i] = (edge_lens[i] / parameters.rest_edge_len) - 1.0);
        let edge_forces = calc_edge_forces(&edge_strains, uevs, parameters.stiffness_edge);
        let avg_tens_strain = edge_strains
            .iter()
            .map(|&es| if es < 0.0 { 0.0 } else { es })
            .sum::<f32>()
            / NVERTS as f32;
        let mut sum_fs = [P2D::default(); NVERTS];
        (0..NVERTS).for_each(|i| {
            sum_fs[i] = rgtp_forces[i] + cyto_forces[i] + edge_forces[i]
                - edge_forces[circ_ix_minus(i as usize, NVERTS)];
        });
        MechState {
            edge_strains,
            rgtp_forces,
            cyto_forces,
            edge_forces,
            avg_tens_strain,
            sum_fs,
        }
    }

    pub fn calc_chem_state(
        &self,
        geom_state: &GeomState,
        mech_state: &MechState,
        rac_rand_state: &RacRandState,
        inter_state: &InteractionState,
        parameters: &Parameters,
    ) -> ChemState {
        let GeomState { edge_lens, .. } = geom_state;
        let mut avg_edge_lens: [f32; NVERTS] = [0.0_f32; NVERTS];
        (0..NVERTS).for_each(|i| {
            let im1 = circ_ix_minus(i as usize, NVERTS);
            avg_edge_lens[i] = (edge_lens[i] + edge_lens[im1]) / 2.0;
        });

        let conc_rac_acts = calc_conc_rgtps(&avg_edge_lens, &self.rac_acts);
        let conc_rac_inacts = calc_conc_rgtps(&avg_edge_lens, &self.rac_inacts);
        let conc_rho_acts = calc_conc_rgtps(&avg_edge_lens, &self.rho_acts);
        let conc_rho_inacts = calc_conc_rgtps(&avg_edge_lens, &self.rho_inacts);

        let kgtps_rac = calc_kgtps_rac(
            &self.rac_acts,
            &conc_rac_acts,
            &rac_rand_state.x_rands,
            &inter_state.x_coas,
            &inter_state.x_chemoas,
            parameters.kgtp_rac,
            parameters.kgtp_rac_auto,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let MechState {
            avg_tens_strain, ..
        } = mech_state;
        let x_tens = parameters.tension_inhib
            * hill_function3(parameters.halfmax_tension_inhib, *avg_tens_strain);
        let kdgtps_rac = calc_kdgtps_rac(
            &self.rac_acts,
            &conc_rho_acts,
            &inter_state.x_cils,
            x_tens,
            parameters.kdgtp_rac,
            parameters.kdgtp_rho_on_rac,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let kgtps_rho = calc_kgtps_rho(
            &self.rho_acts,
            &conc_rho_acts,
            &inter_state.x_cils,
            parameters.kgtp_rho,
            parameters.halfmax_vertex_rgtp_conc,
            parameters.kgtp_rho_auto,
        );
        let kdgtps_rho = calc_kdgtps_rho(
            &self.rho_acts,
            &conc_rac_acts,
            parameters.kdgtp_rho,
            parameters.kdgtp_rac_on_rho,
            parameters.halfmax_vertex_rgtp_conc,
        );
        let rac_act_net_fluxes =
            calc_net_fluxes(&edge_lens, parameters.diffusion_rgtp, &conc_rac_acts);
        let rho_act_net_fluxes =
            calc_net_fluxes(&edge_lens, parameters.diffusion_rgtp, &conc_rho_acts);
        let rac_inact_net_fluxes =
            calc_net_fluxes(&edge_lens, parameters.diffusion_rgtp, &conc_rac_inacts);
        let rho_inact_net_fluxes =
            calc_net_fluxes(&edge_lens, parameters.diffusion_rgtp, &conc_rho_inacts);

        let rac_cyto = parameters.total_rgtp
            - self.rac_acts.iter().sum::<f32>()
            - self.rac_inacts.iter().sum::<f32>();
        let rho_cyto = parameters.total_rgtp
            - self.rho_acts.iter().sum::<f32>()
            - self.rho_inacts.iter().sum::<f32>();
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

    pub fn dynamics_f(
        state: &CoreState,
        rac_rand_state: &RacRandState,
        inter_state: &InteractionState,
        world_parameters: &WorldParameters,
        parameters: &Parameters,
    ) -> CoreState {
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
        let mut delta = CoreState::default();
        for i in 0..NVERTS {
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
            delta.vertex_coords[i] = (1.0 / world_parameters.vertex_eta) * mech_state.sum_fs[i];
        }
        delta
    }

    pub fn new(
        vertex_coords: [P2D; NVERTS],
        init_rac: RgtpDistribution,
        init_rho: RgtpDistribution,
    ) -> CoreState {
        // x_cils: [f32; NVERTS], x_coas: [f32; NVERTS], x_chemoas: [f32; NVERTS], x_rands: [f32; NVERTS], x_bdrys: [f32; NVERTS];
        CoreState {
            vertex_coords,
            rac_acts: init_rac.active,
            rac_inacts: init_rac.inactive,
            rho_acts: init_rho.active,
            rho_inacts: init_rho.inactive,
        }
    }

    pub fn scalar_mul(&self, s: f32) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = s * self.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] * s;
            rac_inacts[i] = self.rac_inacts[i] * s;
            rho_acts[i] = self.rho_acts[i] * s;
            rho_inacts[i] = self.rho_inacts[i] * s;
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn scalar_add(&self, s: f32) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = s + self.vertex_coords[i];
            rac_acts[i] = self.rac_acts[i] + s;
            rac_inacts[i] = self.rac_inacts[i] + s;
            rho_acts[i] = self.rho_acts[i] + s;
            rho_inacts[i] = self.rho_inacts[i] + s;
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn abs(&self) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].abs();
            rac_acts[i] = self.rac_acts[i].abs();
            rac_inacts[i] = self.rac_inacts[i].abs();
            rho_acts[i] = self.rho_acts[i].abs();
            rho_inacts[i] = self.rho_inacts[i].abs();
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn powi(&self, x: i32) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].powi(x);
            rac_acts[i] = self.rac_acts[i].powi(x);
            rac_inacts[i] = self.rac_inacts[i].powi(x);
            rho_acts[i] = self.rho_acts[i].powi(x);
            rho_inacts[i] = self.rho_inacts[i].powi(x);
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn max(&self, other: &CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].max(&other.vertex_coords[i]);
            rac_acts[i] = max_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] = max_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] = max_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] = max_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn min(&self, other: &CoreState) -> CoreState {
        let mut vertex_coords = [P2D::default(); NVERTS];
        let mut rac_acts = [0.0_f32; NVERTS];
        let mut rac_inacts = [0.0_f32; NVERTS];
        let mut rho_acts = [0.0_f32; NVERTS];
        let mut rho_inacts = [0.0_f32; NVERTS];

        for i in 0..(NVERTS) {
            vertex_coords[i] = vertex_coords[i].min(&other.vertex_coords[i]);
            rac_acts[i] = min_f32(self.rac_acts[i], other.rac_acts[i]);
            rac_inacts[i] = min_f32(self.rac_inacts[i], other.rac_inacts[i]);
            rho_acts[i] = min_f32(self.rho_acts[i], other.rho_acts[i]);
            rho_inacts[i] = min_f32(self.rho_inacts[i], other.rho_inacts[i]);
        }

        CoreState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
        }
    }

    pub fn sum(&self) -> f32 {
        let mut r: f32 = 0.0;

        for i in 0..(NVERTS) {
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

    pub fn validate(&self, loc_str: &str, parameters: &Parameters) {
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
        if sum_rac_mem > parameters.total_rgtp || sum_rac_mem < 0.0 {
            panic!("{}: problem in sum of rac_mem: {}", loc_str, sum_rac_mem);
        }
        let sum_rho_mem = self.rho_inacts.iter().sum::<f32>() + self.rho_acts.iter().sum::<f32>();
        if sum_rho_mem > parameters.total_rgtp || sum_rho_mem < 0.0 {
            panic!("{}: problem in sum of rho_mem: {}", loc_str, sum_rho_mem);
        }
        // println!("{}: successfully validated", loc_str)
    }
}
