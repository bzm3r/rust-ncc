// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::chemistry::{
    calc_conc_rgtps, calc_k_mem_offs, calc_k_mem_on, calc_kdgtps_rac, calc_kdgtps_rho,
    calc_kgtps_rac, calc_kgtps_rho, calc_net_fluxes, gen_rgtp_distrib, RgtpLayout,
};
use crate::consts::NVERTS;
use crate::math::{calc_poly_area, hill_function, Radians, P2D};
use crate::mechanics::{
    calc_cyto_forces, calc_edge_forces, calc_edge_lens, calc_edge_strains,
    calc_edge_unit_vecs, calc_global_strain, calc_rgtp_forces,
};
use crate::parameters::Parameters;
use crate::utils::{circ_ix_minus, circ_ix_plus};
use std::f32::consts::PI;

#[derive(Clone, Debug)]
pub struct CellState {
    vertex_coords: Vec<P2D>,
    rac_acts: Vec<f32>,
    rac_inacts: Vec<f32>,
    rho_acts: Vec<f32>,
    rho_inacts: Vec<f32>,
    x_cils: Vec<f32>,
    x_chemoas: Vec<f32>,
    x_coas: Vec<f32>,
    x_rands: Vec<f32>,
    x_bdrys: Vec<f32>,
}

pub struct DependentVars {
    pub kdgtps_rac: Vec<f32>,
    pub kgtps_rac: Vec<f32>,
    pub rac_act_net_fluxes: Vec<f32>,
    pub rac_inact_net_fluxes: Vec<f32>,
    pub kdgtps_rho: Vec<f32>,
    pub kgtps_rho: Vec<f32>,
    pub rho_act_net_fluxes: Vec<f32>,
    pub rho_inact_net_fluxes: Vec<f32>,
    pub rgtp_forces: Vec<P2D>,
    pub cyto_forces: Vec<P2D>,
    pub edge_forces: Vec<P2D>,
    pub rac_mem_on: f32,
    pub rac_mem_offs: Vec<f32>,
    pub rho_mem_on: f32,
    pub rho_mem_offs: Vec<f32>,
}

fn increment_f32s(xs: &[f32], dxs: &[f32]) -> Vec<f32> {
    xs.iter().zip(dxs.iter()).map(|(x, dx)| x + dx).collect()
}

fn increment_vec2ds(xs: &[P2D], dxs: &[P2D]) -> Vec<P2D> {
    xs.iter().zip(dxs.iter()).map(|(x, dx)| x + dx).collect()
}

impl CellState {
    fn calc_dep_vars(&self, params: &Parameters) -> DependentVars {
        let edge_uvs = calc_edge_unit_vecs(&self.vertex_coords);
        let edge_lens: Vec<f32> = calc_edge_lens(&edge_uvs);
        let avg_edge_lens: Vec<f32> = (0..NVERTS)
            .map(|i| {
                let im1 = circ_ix_minus(i, NVERTS);
                (edge_lens[i] + edge_lens[im1]) / 2.0
            })
            .collect();
        let unit_inward_vecs: Vec<P2D> = (0..NVERTS)
            .map(|i| {
                let im1 = circ_ix_minus(i, NVERTS);
                let tangent = (edge_uvs[i] + edge_uvs[im1]).unitize();
                tangent.normal()
            })
            .collect();
        let conc_rac_acts = calc_conc_rgtps(&avg_edge_lens, &self.rac_acts);
        let conc_rac_inacts = calc_conc_rgtps(&avg_edge_lens, &self.rac_inacts);
        let conc_rho_acts = calc_conc_rgtps(&avg_edge_lens, &self.rho_acts);
        let conc_rho_inacts = calc_conc_rgtps(&avg_edge_lens, &self.rho_inacts);

        let kgtps_rac = calc_kgtps_rac(
            &self.rac_acts,
            &conc_rac_acts,
            &self.x_rands,
            &self.x_coas,
            &self.x_chemoas,
            params.kgtp_rac,
            params.kgtp_rac_auto,
            params.halfmax_vertex_rgtp_conc,
        );
        let x_tens = hill_function(
            params.halfmax_tension_inhib,
            calc_global_strain(&edge_lens, params.rest_edge_len, NVERTS),
        );
        let kdgtps_rac = calc_kdgtps_rac(
            &self.rac_acts,
            &conc_rac_acts,
            &self.x_cils,
            x_tens,
            params.kdgtp_rac,
            params.kdgtp_rho_on_rac,
            params.halfmax_vertex_rgtp_conc,
        );
        let kgtps_rho = calc_kgtps_rho(
            &self.rho_acts,
            &conc_rac_acts,
            &self.x_cils,
            params.kgtp_rho,
            params.halfmax_vertex_rgtp_conc,
            params.kgtp_rho_auto,
        );
        let kdgtps_rho = calc_kdgtps_rho(
            &self.rho_acts,
            &conc_rac_acts,
            params.kdgtp_rho,
            params.kdgtp_rac_on_rho,
            params.halfmax_vertex_rgtp_conc,
        );
        let rac_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            &avg_edge_lens,
            params.diffusion_rgtp,
            &conc_rac_acts,
        );
        let rho_act_net_fluxes = calc_net_fluxes(
            &edge_lens,
            &avg_edge_lens,
            params.diffusion_rgtp,
            &conc_rho_acts,
        );
        let rac_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            &avg_edge_lens,
            params.diffusion_rgtp,
            &conc_rac_inacts,
        );
        let rho_inact_net_fluxes = calc_net_fluxes(
            &edge_lens,
            &avg_edge_lens,
            params.diffusion_rgtp,
            &conc_rho_inacts,
        );

        let rac_cyto =
            1.0_f32 - self.rac_acts.iter().sum::<f32>() - self.rac_inacts.iter().sum::<f32>();
        let rac_mem_on = calc_k_mem_on(rac_cyto, params.k_mem_on);
        let rac_mem_offs = calc_k_mem_offs(&self.rho_inacts, params.k_mem_off);
        let rho_cyto =
            1.0 - self.rho_acts.iter().sum::<f32>() - self.rho_inacts.iter().sum::<f32>();
        let rho_mem_on = calc_k_mem_on(rho_cyto, params.k_mem_on);
        let rho_mem_offs = calc_k_mem_offs(&self.rac_inacts, params.k_mem_off);
        let rgtp_forces = calc_rgtp_forces(
            &self.rac_acts,
            &self.rho_acts,
            &unit_inward_vecs,
            params.max_protrusive_f,
            params.max_retractive_f,
            params.vertex_rgtp_act_at_max_f,
        );
        let cyto_forces = calc_cyto_forces(
            &self.vertex_coords,
            &unit_inward_vecs,
            params.rest_area,
            params.stiffness_ctyo,
        );
        let edge_strains = calc_edge_strains(&edge_lens, params.rest_edge_len);
        let edge_forces = calc_edge_forces(&edge_strains, &edge_uvs, params.stiffness_edge);

        DependentVars {
            kdgtps_rac,
            kgtps_rac,
            rac_act_net_fluxes,
            rac_inact_net_fluxes,
            kdgtps_rho,
            kgtps_rho,
            rho_act_net_fluxes,
            rho_inact_net_fluxes,
            rac_mem_on,
            rac_mem_offs,
            rho_mem_on,
            rho_mem_offs,
            rgtp_forces,
            cyto_forces,
            edge_forces,
        }
    }

    fn euler_integrate(&self, params: &Parameters, dt: f32) -> CellState {
        let DependentVars {
            kdgtps_rac,
            kgtps_rac,
            rac_act_net_fluxes,
            rac_inact_net_fluxes,
            kdgtps_rho,
            kgtps_rho,
            rho_act_net_fluxes,
            rho_inact_net_fluxes,
            rac_mem_on,
            rac_mem_offs,
            rho_mem_on,
            rho_mem_offs,
            rgtp_forces,
            cyto_forces,
            edge_forces,
        } = self.calc_dep_vars(params);
        let mut delta_rac_acts = vec![0.0; NVERTS];
        let mut delta_rac_inacts = vec![0.0; NVERTS];
        let mut delta_rho_acts = vec![0.0; NVERTS];
        let mut delta_rho_inacts = vec![0.0; NVERTS];
        let mut delta_pos = vec![P2D::default(); NVERTS];
        for i in 0..NVERTS {
            let inactivated_rac = kdgtps_rac[i] * self.rac_acts[i];
            let activated_rac = kgtps_rac[i] * self.rac_inacts[i];
            let delta_rac_activated = activated_rac - inactivated_rac;
            delta_rac_acts[i] = (delta_rac_activated + rac_act_net_fluxes[i]) * dt;
            delta_rac_inacts[i] =
                (-1.0 * delta_rac_activated + rac_inact_net_fluxes[i] + rac_mem_on - rac_mem_offs[i]) * dt;

            let inactivated_rho = kdgtps_rho[i] * self.rho_acts[i];
            let activated_rho = kgtps_rho[i] * self.rho_inacts[i];
            let delta_rho_activated = activated_rho - inactivated_rho;
            delta_rho_acts[i] = (delta_rho_activated + rho_act_net_fluxes[i]) * dt;
            delta_rho_inacts[i] =
                (-1.0 * delta_rho_activated + rho_inact_net_fluxes[i] + rho_mem_on - rho_mem_offs[i]) * dt;

            delta_pos[i] = ((rgtp_forces[i] + cyto_forces[i] + edge_forces[i]
                - edge_forces[circ_ix_minus(i, NVERTS)])
                / params.vertex_eta) * dt;
        }

        let vertex_coords = increment_vec2ds(&self.vertex_coords, &delta_pos);

        let rac_acts = increment_f32s(&self.rac_acts, &delta_rac_acts);
        let rac_inacts = increment_f32s(&self.rac_inacts, &delta_rac_inacts);
        let rho_acts = increment_f32s(&self.rho_acts, &delta_rho_acts);
        let rho_inacts = increment_f32s(&self.rho_inacts, &delta_rho_inacts);

        CellState::new(vertex_coords, rac_acts, rac_inacts, rho_acts, rho_inacts)
    }

    pub fn new(
        vertex_coords: Vec<P2D>,
        rac_acts: Vec<f32>,
        rac_inacts: Vec<f32>,
        rho_acts: Vec<f32>,
        rho_inacts: Vec<f32>,
    ) -> CellState {
        // x_cils: Vec<f32>, x_coas: Vec<f32>, x_chemoas: Vec<f32>, x_rands: Vec<f32>, x_bdrys: Vec<f32>
        CellState {
            vertex_coords,
            rac_acts,
            rac_inacts,
            rho_acts,
            rho_inacts,
            x_cils: vec![0.0; NVERTS],
            x_chemoas: vec![0.0; NVERTS],
            x_coas: vec![0.0; NVERTS],
            x_rands: vec![0.0; NVERTS],
            x_bdrys: vec![0.0; NVERTS],
        }
    }
}

/// Model cell.
pub struct Cell {
    /// Cell index.
    pub(crate) ix: usize,
    /// Index of cell type.
    // group_ix: usize,
    // /// If `true`, `evolve` does not change current cell state.
    // skip_dynamics: bool,
    /// Parameters.
    params: Parameters,
    /// Current state.
    pub(crate) state: CellState,
    /// Last state.
    history: Option<CellState>,
    pub group_ix: usize,
}

const NCELLS: usize = 16;

fn gen_vertex_coords(centroid: P2D, radius: f32) -> Vec<P2D> {
    (0..NVERTS).map(|vix| {
        let vf = (vix as f32) / (NVERTS as f32);
        let theta = 2.0 * PI * vf;
        P2D {
            x: &centroid.x + theta.cos() * radius,
            y: &centroid.y + theta.sin() * radius,
        }
    }).collect()
}

impl Cell {
    pub fn new(ix: usize, group_ix: usize, params: Parameters, centroid: P2D) -> Cell {
        let vertex_coords = gen_vertex_coords(centroid, params.cell_r);
        let (rac_acts, rac_inacts) = gen_rgtp_distrib(
            &vertex_coords,
            params.init_frac_active,
            params.init_frac_inactive,
            RgtpLayout::Random,
        );

        let (rho_acts, rho_inacts) = gen_rgtp_distrib(
            &vertex_coords,
            params.init_frac_active,
            params.init_frac_inactive,
            RgtpLayout::Random,
        );

        Cell {
            ix,
            group_ix,
            params,
            state: CellState::new(vertex_coords, rac_acts, rac_inacts, rho_acts, rho_inacts),
            history: None,
        }
    }

    pub fn simulate(&mut self) {
        self.history = Some(self.state.clone());
        let nsteps: u32 = 10;
        let dt = 1.0 / (nsteps as f32);
        for _ in 0..nsteps {
            self.state = self.state.euler_integrate(&self.params, dt);
        }
    }
}
