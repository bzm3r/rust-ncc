// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
use crate::math::{P2D, Radians, hill_function};
use crate::consts::NVERTS;
use crate::utils::circ_ix_plus;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use crate::parameters::Parameters;
use crate::mechanics::calc_global_strain;


pub struct BiasDefn {
    range: (Radians, Radians),
    //strength: f32,
}

pub enum RgtpLayout {
    Random,
    //Biased(BiasDefn),
}

impl RgtpLayout {
    pub fn gen_distrib(&self, vertex_coords: &[P2D]) -> [f32; NVERTS] {
        let mut rgtp: [f32; NVERTS] = [0.0; NVERTS];
        match self {
            RgtpLayout::Random => {
                let distr: Uniform<f32> = Uniform::new_inclusive(0.0, 1.0);
                let mut rng = thread_rng();
                rgtp.iter_mut().for_each(|e| {
                    *e = rng.sample(distr);
                });
            }
            // RgtpLayout::Biased(bd) => {
            //     rgtp.iter_mut().enumerate().for_each(|(i, e)| {
            //         if vertex_coords[i].direction().between(bd.range.0, bd.range.1) {
            //             *e = 1.0;
            //         }
            //     });
            // }
        }
        let sum: f32 = rgtp.iter().sum();
        rgtp.iter_mut().for_each(|e| *e = *e / sum);
        rgtp
    }
}

pub fn gen_rgtp_distrib(
    vertex_coords: &[P2D],
    active_frac: f32,
    inactive_frac: f32,
    bias: RgtpLayout,
) -> (Vec<f32>, Vec<f32>) {
    let distrib = bias.gen_distrib(vertex_coords);
    let n = distrib.len();
    let active_frac = active_frac / (n as f32);
    let inactive_frac = inactive_frac / (n as f32);

    let active = distrib
        .iter()
        .map(|&x| active_frac * x)
        .collect::<Vec<f32>>();
    let inactive = distrib
        .iter()
        .map(|&x| inactive_frac * x)
        .collect::<Vec<f32>>();

    (active, inactive)
}

fn calc_directed_fluxes(edge_lens: &[f32], avg_edge_lens: &[f32], rgtp_d: f32, conc_rgtps: &[f32]) -> Vec<f32> {
    let nvs = conc_rgtps.len();
    conc_rgtps
        .iter()
        .enumerate()
        .map(|(i, &src_conc)| {
            let plus_i = circ_ix_plus(i, nvs);
            let dst_conc = conc_rgtps[plus_i];
            let dst_avg_el = avg_edge_lens[plus_i];
            let src_avg_el = avg_edge_lens[i];
            let el_src_to_dst = edge_lens[i];
            rgtp_d * ((dst_conc / dst_avg_el) + (src_conc / src_avg_el)) / el_src_to_dst
        })
        .collect()
}

pub fn calc_net_fluxes(edge_lens: &[f32], avg_edge_lens: &[f32], rgtp_d: f32, conc_rgtps: &[f32]) -> Vec<f32> {
    let directed_fluxes = calc_directed_fluxes(edge_lens, avg_edge_lens, rgtp_d, conc_rgtps);
    let nvs = directed_fluxes.len();
    (0..directed_fluxes.len())
        .map(|i| {
            let min_i = circ_ix_plus(i, nvs);
            directed_fluxes[min_i] - directed_fluxes[i]
        })
        .collect()
}


pub fn calc_conc_rgtps(avg_edge_lens: &[f32], rgtps: &[f32]) -> Vec<f32> {
    rgtps
        .iter()
        .zip(avg_edge_lens.iter())
        .map(|(&r, &ael)| r / ael)
        .collect()
}

pub fn calc_kgtps_rac(
    rac_acts: &[f32],
    conc_rac_acts: &[f32],
    x_rands: &[f32],
    x_coas: &[f32],
    x_chemoas: &[f32],
    kgtp_rac_base: f32,
    kgtp_rac_auto: f32,
    halfmax_rac_conc: f32
) -> Vec<f32> {
    let nvs = rac_acts.len();
    let mut kgtps_rac = vec![0.0_f32; nvs];

    for i in 0..nvs {
        let base = (x_rands[i] + x_coas[i] + 1.0) * kgtp_rac_base;
        let auto = hill_function(halfmax_rac_conc, conc_rac_acts[i])
            * (1.0 + x_chemoas[i])
            * kgtp_rac_auto;
        kgtps_rac[i] = base + auto;
    }

    kgtps_rac
}

pub fn calc_kdgtps_rac(
    rac_acts: &[f32],
    conc_rho_acts: &[f32],
    x_cils: &[f32],
    x_tens: f32,
    kdgtp_rac_base: f32,
    kdgtp_rho_on_rac: f32,
    halfmax_conc_rho: f32,
) -> Vec<f32> {
    let nvs = rac_acts.len();
    let mut kdgtps_rac = vec![0.0_f32; nvs];

    for i in 0..nvs {
        let base = (1.0 + x_tens + x_cils[i]) * kdgtp_rac_base;
        let mutual = hill_function(halfmax_conc_rho, conc_rho_acts[i]) * kdgtp_rho_on_rac;
        kdgtps_rac[i] = base + mutual;
    }

    kdgtps_rac
}

pub fn calc_kgtps_rho(
    rho_acts: &[f32],
    conc_rho_acts: &[f32],
    x_cils: &[f32],
    kgtp_rho_base: f32,
    halfmax_rho_thresh: f32,
    kgtp_rho_auto: f32,
) -> Vec<f32> {
    let nvs = rho_acts.len();
    let mut kgtps_rho = vec![0.0_f32; nvs];

    for i in 0..nvs {
        let base = (1.0 + x_cils[i]) * kgtp_rho_base;
        let auto = hill_function(halfmax_rho_thresh, conc_rho_acts[i]) * kgtp_rho_auto;
        kgtps_rho[i] = base + auto;
    }

    kgtps_rho
}

pub fn calc_kdgtps_rho(
    rho_acts: &[f32],
    conc_rac_acts: &[f32],
    kdgtp_rho_base: f32,
    kdgtp_rac_on_rho: f32,
    halfmax_conc_rac: f32,
) -> Vec<f32> {
    let nvs = rho_acts.len();
    let mut kdgtps_rho = vec![0.0_f32; nvs];

    for i in 0..nvs {
        let mutual = hill_function(halfmax_conc_rac, conc_rac_acts[i]) * kdgtp_rac_on_rho;
        kdgtps_rho[i] = kdgtp_rho_base + mutual;
    }

    kdgtps_rho
}

pub fn calc_k_mem_on(cytosol_frac: f32, k_mem_on: f32) -> f32 {
    (cytosol_frac * k_mem_on) / (NVERTS as f32)
}

pub fn calc_k_mem_offs(inacts: &[f32], k_mem_off: f32) -> Vec<f32> {
    inacts.iter().map(|&inact| k_mem_off * inact).collect()
}
