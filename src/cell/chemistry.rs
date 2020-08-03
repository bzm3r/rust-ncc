// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::hill_function3;
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;
use avro_schema_derive::Schematize;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};

pub enum RgtpLayout {
    Random,
    BiasedVertices(Vec<usize>),
    //Biased(BiasDefn),
}

impl RgtpLayout {
    pub fn gen_distrib(&self, frac_to_distrib: f32) -> [f32; NVERTS as usize] {
        let mut rgtp_distrib: [f32; NVERTS as usize] = [0.0_f32; NVERTS as usize];
        match self {
            RgtpLayout::Random => {
                let distrib: Uniform<f32> = Uniform::new_inclusive(0.0, 1.0);
                let mut rng = thread_rng();
                rgtp_distrib.iter_mut().for_each(|e| {
                    *e = rng.sample(distrib);
                });
            }
            RgtpLayout::BiasedVertices(verts) => {
                (0..NVERTS as usize).for_each(|ix| {
                    if let Some(_) = verts.iter().find(|&&v| ix == v) {
                        rgtp_distrib[ix] = 1.0;
                    }
                });
            }
        }
        let sum: f32 = rgtp_distrib.iter().sum();
        rgtp_distrib
            .iter_mut()
            .for_each(|e| *e = *e * frac_to_distrib / sum);
        rgtp_distrib
    }
}

pub fn gen_rgtp_distrib(
    active_frac: f32,
    inactive_frac: f32,
    rgtp_layout: &RgtpLayout,
) -> ([f32; NVERTS as usize], [f32; NVERTS as usize]) {
    (
        rgtp_layout.gen_distrib(active_frac),
        rgtp_layout.gen_distrib(inactive_frac),
    )
}

fn calc_directed_fluxes(
    edge_lens: &[f32; NVERTS as usize],
    rgtp_d: f32,
    conc_rgtps: &[f32; NVERTS as usize],
) -> [f32; NVERTS as usize] {
    let mut r = [0.0_f32; NVERTS as usize];
    for i in 0..NVERTS as usize {
        let plus_i = circ_ix_plus(i, NVERTS as usize);
        r[i] = -1.0 * rgtp_d * (conc_rgtps[plus_i] - conc_rgtps[i]) / edge_lens[i];
    }
    r
}

pub fn calc_net_fluxes(
    edge_lens: &[f32; NVERTS as usize],
    rgtp_d: f32,
    conc_rgtps: &[f32; NVERTS as usize],
) -> [f32; NVERTS as usize] {
    let directed_fluxes = calc_directed_fluxes(edge_lens, rgtp_d, conc_rgtps);
    let mut r = [0.0_f32; NVERTS as usize];
    (0..NVERTS as usize).for_each(|i| {
        let min_i = circ_ix_minus(i, NVERTS as usize);
        r[i] = directed_fluxes[min_i] - directed_fluxes[i];
    });
    r
}

pub fn calc_conc_rgtps(
    avg_edge_lens: &[f32; NVERTS as usize],
    rgtps: &[f32; NVERTS as usize],
) -> [f32; NVERTS as usize] {
    let mut r = [0.0_f32; NVERTS as usize];
    (0..NVERTS as usize).for_each(|i| r[i] = rgtps[i] / avg_edge_lens[i]);
    r
}

#[allow(clippy::too_many_arguments)]
pub fn calc_kgtps_rac(
    rac_acts: &[f32; NVERTS as usize],
    conc_rac_acts: &[f32; NVERTS as usize],
    x_rands: &[f32; NVERTS as usize],
    x_coas: &[f32; NVERTS as usize],
    x_chemoas: &[f32; NVERTS as usize],
    kgtp_rac_base: f32,
    kgtp_rac_auto: f32,
    halfmax_rac_conc: f32,
) -> [f32; NVERTS as usize] {
    let nvs = rac_acts.len();
    let mut kgtps_rac = [0.0_f32; NVERTS as usize];

    for i in 0..nvs {
        let base = (x_rands[i] + x_coas[i] + 1.0) * kgtp_rac_base;
        let auto = hill_function3(halfmax_rac_conc, conc_rac_acts[i])
            * (1.0 + x_chemoas[i])
            * kgtp_rac_auto;
        kgtps_rac[i] = base + auto;
    }

    kgtps_rac
}

pub fn calc_kdgtps_rac(
    rac_acts: &[f32; NVERTS as usize],
    conc_rho_acts: &[f32; NVERTS as usize],
    x_cils: &[f32; NVERTS as usize],
    x_tens: f32,
    kdgtp_rac_base: f32,
    kdgtp_rho_on_rac: f32,
    halfmax_conc_rho: f32,
) -> [f32; NVERTS as usize] {
    let nvs = rac_acts.len();
    let mut kdgtps_rac = [0.0_f32; NVERTS as usize];

    for i in 0..nvs {
        let base = (1.0 + x_tens + x_cils[i]) * kdgtp_rac_base;
        let mutual = hill_function3(halfmax_conc_rho, conc_rho_acts[i]) * kdgtp_rho_on_rac;
        kdgtps_rac[i] = base + mutual;
    }

    kdgtps_rac
}

pub fn calc_kgtps_rho(
    rho_acts: &[f32; NVERTS as usize],
    conc_rho_acts: &[f32; NVERTS as usize],
    x_cils: &[f32; NVERTS as usize],
    kgtp_rho_base: f32,
    halfmax_rho_thresh: f32,
    kgtp_rho_auto: f32,
) -> [f32; NVERTS as usize] {
    let nvs = rho_acts.len();
    let mut kgtps_rho = [0.0_f32; NVERTS as usize];

    for i in 0..nvs {
        let base = (1.0 + x_cils[i]) * kgtp_rho_base;
        let auto = hill_function3(halfmax_rho_thresh, conc_rho_acts[i]) * kgtp_rho_auto;
        kgtps_rho[i] = base + auto;
    }

    kgtps_rho
}

pub fn calc_kdgtps_rho(
    rho_acts: &[f32; NVERTS as usize],
    conc_rac_acts: &[f32; NVERTS as usize],
    kdgtp_rho_base: f32,
    kdgtp_rac_on_rho: f32,
    halfmax_conc_rac: f32,
) -> [f32; NVERTS as usize] {
    let nvs = rho_acts.len();
    let mut kdgtps_rho = [0.0_f32; NVERTS as usize];

    for i in 0..nvs {
        let mutual = hill_function3(halfmax_conc_rac, conc_rac_acts[i]) * kdgtp_rac_on_rho;
        kdgtps_rho[i] = kdgtp_rho_base + mutual;
    }

    kdgtps_rho
}

#[derive(Clone, Deserialize, Serialize, Schematize)]
pub struct RacRandomState {
    pub x_rands: [f32; NVERTS as usize],
}

impl RacRandomState {
    pub fn init() -> RacRandomState {
        RacRandomState {
            x_rands: [0.0; NVERTS as usize],
        }
    }

    pub fn update(&self, _tstep: u32) -> RacRandomState {
        Self::init()
    }
}
