// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::interactions::RgtpState;
use crate::math::hill_function3;
use crate::model_cell::core_state::fmt_var_arr;
use crate::parameters::Parameters;
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::world::RandomEventGenerator;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Display;

pub enum DistributionType {
    Random,
    Specific([bool; NVERTS as usize]),
}

pub struct DistributionScheme {
    pub frac: f32,
    pub ty: DistributionType,
}

impl DistributionScheme {
    fn scaled_unitize(frac: f32, mut distrib: [f32; NVERTS as usize]) -> [f32; NVERTS as usize] {
        let sum: f32 = distrib.iter().sum();
        distrib.iter_mut().for_each(|e| *e = *e * frac / sum);
        distrib
    }

    fn gen_random(frac: f32) -> [f32; NVERTS as usize] {
        let mut r = [0.0; NVERTS as usize];
        let distrib: Uniform<f32> = Uniform::new_inclusive(0.0, 1.0);
        let mut rng = thread_rng();
        r.iter_mut().for_each(|e| {
            *e = rng.sample(distrib);
        });
        Self::scaled_unitize(frac, r)
    }

    fn gen_specific(frac: f32, verts: &[bool; NVERTS as usize]) -> [f32; NVERTS as usize] {
        let mut r = [0.0; NVERTS as usize];
        verts.iter().zip(r.iter_mut()).for_each(|(&marked, e)| {
            if marked {
                *e = 1.0;
            }
        });
        Self::scaled_unitize(frac, r)
    }

    pub fn generate(&self) -> [f32; NVERTS] {
        match &self.ty {
            DistributionType::Random => Self::gen_random(self.frac),
            DistributionType::Specific(marks) => Self::gen_specific(self.frac, marks),
        }
    }
}

#[allow(unused)]
#[derive(Clone, Copy)]
pub struct RgtpDistribution {
    pub active: [f32; NVERTS as usize],
    pub inactive: [f32; NVERTS as usize],
}

impl RgtpDistribution {
    pub fn generate(
        act_distrib: DistributionScheme,
        inact_distrib: DistributionScheme,
    ) -> Result<RgtpDistribution, String> {
        if act_distrib.frac + inact_distrib.frac > 1.0 {
            Err("active + inactive > 1.0".to_string())
        } else {
            Ok(RgtpDistribution {
                active: act_distrib.generate(),
                inactive: inact_distrib.generate(),
            })
        }
    }
}

fn calc_directed_fluxes(
    edge_lens: &[f32; NVERTS],
    rgtp_d: f32,
    conc_rgtps: &[f32; NVERTS],
) -> [f32; NVERTS] {
    let mut r = [0.0_f32; NVERTS];
    for i in 0..NVERTS {
        let plus_i = circ_ix_plus(i, NVERTS);
        r[i] = -1.0 * rgtp_d * (conc_rgtps[plus_i] - conc_rgtps[i]) / edge_lens[i];
    }
    r
}

pub fn calc_net_fluxes(
    edge_lens: &[f32; NVERTS],
    rgtp_d: f32,
    conc_rgtps: &[f32; NVERTS],
) -> [f32; NVERTS] {
    let directed_fluxes = calc_directed_fluxes(edge_lens, rgtp_d, conc_rgtps);
    let mut r = [0.0_f32; NVERTS];
    (0..NVERTS).for_each(|i| {
        let min_i = circ_ix_minus(i, NVERTS);
        r[i] = directed_fluxes[min_i] - directed_fluxes[i];
    });
    r
}

pub fn calc_conc_rgtps(avg_edge_lens: &[f32; NVERTS], rgtps: &[f32; NVERTS]) -> [f32; NVERTS] {
    let mut r = [0.0_f32; NVERTS];
    (0..NVERTS).for_each(|i| r[i] = rgtps[i] / avg_edge_lens[i]);
    r
}

#[allow(clippy::too_many_arguments)]
pub fn calc_kgtps_rac(
    rac_acts: &[f32; NVERTS],
    conc_rac_acts: &[f32; NVERTS],
    x_rands: &[f32; NVERTS],
    x_coas: &[f32; NVERTS],
    x_chemoas: &[f32; NVERTS],
    x_cals: &[f32; NVERTS],
    kgtp_rac_base: f32,
    kgtp_rac_auto: f32,
    halfmax_rac_conc: f32,
) -> [f32; NVERTS] {
    let nvs = rac_acts.len();
    let mut kgtps_rac = [0.0_f32; NVERTS];

    for i in 0..nvs {
        let base = (x_cals[i] + x_rands[i] + x_coas[i] + 1.0) * kgtp_rac_base;
        let auto_factor = {
            let af = hill_function3(halfmax_rac_conc, conc_rac_acts[i]) * (1.0 + x_chemoas[i]);
            if af > 1.25 {
                1.25
            } else {
                af
            }
        };
        let auto = auto_factor * kgtp_rac_auto;
        kgtps_rac[i] = base + auto;
    }

    kgtps_rac
}

pub fn calc_kdgtps_rac(
    rac_acts: &[f32; NVERTS],
    conc_rho_acts: &[f32; NVERTS],
    x_cils: &[f32; NVERTS],
    x_tens: f32,
    kdgtp_rac_base: f32,
    kdgtp_rho_on_rac: f32,
    halfmax_conc_rho: f32,
) -> [f32; NVERTS] {
    let nvs = rac_acts.len();
    let mut kdgtps_rac = [0.0_f32; NVERTS];

    for i in 0..nvs {
        // let x_cil =
        //     (x_cils[circ_ix_minus(i, NVERTS)] + x_cils[i] + x_cils[circ_ix_plus(i, NVERTS)]) / 3.0;
        let base = (1.0 + x_tens + x_cils[i]) * kdgtp_rac_base;
        let mutual = hill_function3(halfmax_conc_rho, conc_rho_acts[i]) * kdgtp_rho_on_rac;
        kdgtps_rac[i] = base + mutual;
    }

    kdgtps_rac
}

pub fn calc_kgtps_rho(
    rho_acts: &[f32; NVERTS],
    conc_rho_acts: &[f32; NVERTS],
    x_cils: &[f32; NVERTS],
    kgtp_rho_base: f32,
    halfmax_rho_thresh: f32,
    kgtp_rho_auto: f32,
) -> [f32; NVERTS] {
    let nvs = rho_acts.len();
    let mut kgtps_rho = [0.0_f32; NVERTS];

    for i in 0..nvs {
        // let x_cil =
        //     (x_cils[circ_ix_minus(i, NVERTS)] + x_cils[i] + x_cils[circ_ix_plus(i, NVERTS)]) / 3.0;
        let base = (1.0 + x_cils[i]) * kgtp_rho_base;
        let auto = hill_function3(halfmax_rho_thresh, conc_rho_acts[i]) * kgtp_rho_auto;
        kgtps_rho[i] = base + auto;
    }

    kgtps_rho
}

pub fn calc_kdgtps_rho(
    rho_acts: &[f32; NVERTS],
    conc_rac_acts: &[f32; NVERTS],
    x_cals: &[f32; NVERTS],
    kdgtp_rho_base: f32,
    kdgtp_rac_on_rho: f32,
    halfmax_conc_rac: f32,
) -> [f32; NVERTS] {
    let nvs = rho_acts.len();
    let mut kdgtps_rho = [0.0_f32; NVERTS];

    for i in 0..nvs {
        let mutual = hill_function3(halfmax_conc_rac, conc_rac_acts[i]) * kdgtp_rac_on_rho;
        kdgtps_rho[i] = (1.0 + x_cals[i]) * kdgtp_rho_base + mutual;
    }

    kdgtps_rho
}

#[derive(Copy, Clone, Default, Deserialize, Serialize, Schematize)]
pub struct RacRandState {
    pub next_update: u32,
    pub x_rands: [f32; NVERTS],
}

impl RacRandState {
    pub fn gen_rand_factors(
        rng: &mut SmallRng,
        num_rand_verts: usize,
        rand_mag: f32,
    ) -> [f32; NVERTS] {
        let vs = (0..NVERTS).collect::<Vec<usize>>();
        let mut r = [0.0; NVERTS];
        vs.choose_multiple(rng, num_rand_verts)
            .for_each(|&v| r[v] = rand_mag);
        r
    }

    pub fn init(rng: Option<&mut SmallRng>, parameters: &Parameters) -> RacRandState {
        match rng {
            Some(r) => {
                let ut = Uniform::from(0.0..parameters.rand_avg_t);
                RacRandState {
                    next_update: ut.sample(r).floor() as u32,
                    x_rands: Self::gen_rand_factors(
                        r,
                        parameters.num_rand_vs as usize,
                        parameters.rand_mag,
                    ),
                }
            }
            None => RacRandState {
                next_update: 0,
                x_rands: [0.0; NVERTS],
            },
        }
    }

    pub fn update(
        &self,
        cr: &mut RandomEventGenerator,
        tstep: u32,
        parameters: &Parameters,
    ) -> RacRandState {
        let next_update = tstep + cr.sample().floor() as u32;
        //println!("random update from {} to {}", tstep, next_update);
        let x_rands =
            Self::gen_rand_factors(&mut cr.rng, parameters.num_rand_vs, parameters.rand_mag);
        //println!("{:?}", x_rands);
        RacRandState {
            next_update,
            x_rands,
        }
    }
}

impl Display for RacRandState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_var_arr(f, "rfs", &self.x_rands)?;
        writeln!(f, "next_update: {}", self.next_update)
    }
}

pub fn calc_rgtp_state(
    rac_acts: &[f32; NVERTS],
    rho_acts: &[f32; NVERTS],
    halfmax_rgtp_act: f32,
) -> [RgtpState; NVERTS] {
    let mut r: [RgtpState; NVERTS] = [0.0; NVERTS];
    r.iter_mut()
        .zip(rac_acts.iter().zip(rho_acts.iter()))
        .for_each(|(x, (&rac, &rho))| {
            *x = hill_function3(halfmax_rgtp_act, rac) - hill_function3(halfmax_rgtp_act, rho)
        });
    r
}
