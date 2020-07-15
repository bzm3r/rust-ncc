// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// optension_inhibon. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::consts::NVERTS;
use crate::quantity::{Diffusion, Force, Length, Quantity, Stress, Time, Tinv, Viscosity};
use crate::random::RandomizationType;
use override_derive::Overrides;
use std::f32::consts::PI;

/// Characteristic parameters used for normalization, or elimination of third dimension.
pub struct CharQuants {
    eta: Viscosity,
    f: Force,
    l: Length,
    t: Time,
    l3d: Length,
    k_mem_off: Tinv,
    k_mem_on: Tinv,
    kgtp: Tinv,
    kdgtp: Tinv,
    //coa_dist: Length,
}

impl CharQuants {
    pub fn normalize<T: Quantity>(&self, q: &T) -> f32 {
        let q = q.g();
        let u = q.units();
        (q * self.f.pow(-1.0 * u.f) * self.l.pow(-1.0 * u.l) * self.t.pow(-1.0 * u.t)).value()
    }

    pub fn time(&self) -> f32 {
        self.t.0
    }
}

impl Default for CharQuants {
    fn default() -> Self {
        let eta = {
            let force = Force(290.0_f32).nano().g();
            let time = Time(1.0).g();
            let length = Length(1.0).micro().g();

            (force * time / length)
                .mulf(1.0 / (NVERTS as f32))
                .to_viscosity()
                .unwrap()
        };

        let velocity = Length(3.0).micro().g() / Time(60.0).g();
        let f = (eta.g() * velocity.g()).to_force().unwrap();

        CharQuants {
            eta,
            l: Length(1.0).micro(),
            t: Time(2.0),
            f,
            l3d: Length(10e-6),
            k_mem_off: Tinv(0.02),
            k_mem_on: Tinv(0.15),
            kgtp: Tinv(1e-4),
            kdgtp: Tinv(1e-4),
            //coa_dist: Length(110.0).micro(),
        }
    }
}

#[derive(Overrides)]
pub struct UserParams {
    /// Cell diameter.
    pub cell_diam: Length,
    /// Length criterion to determine if two vertices are "close".
    close_criterion: Length,
    /// Activity at max force.
    rgtp_act_at_max_f: f32,
    /// Stiffness of the membrane-cortex complex.
    stiffness_cortex: Stress,
    /// Typical lamellipod height.
    lm_h: Length,
    /// Halfmax Rho GTPase activity.
    halfmax_rgtp_frac: f32,
    /// Lamellipod stall stress.
    lm_ss: Stress,
    /// Friction force opposing RhoA pulling.
    rho_friction: f32,
    /// Stiffness of cytoplasm.
    stiffness_ctyo: Force,
    /// Diffusion rate of Rho GTPase on membrane as a multiple of characteristic.
    diffusion_rgtp: Diffusion,
    /// Initial fraction of Rho GTPase in cytosol (inactive).
    init_frac_cyto: f32,
    /// Initial fraction of Rho GTPase active.
    init_frac_active: f32,
    /// Hill function exponent for Rho GTPase chemistry.
    h_exp: f32,
    /// Total amount of Rac1 in cell.
    tot_rac: f32,
    /// Total amount of RhoA in cell.
    tot_rho: f32,
    /// Baseline Rac1 activation rate as a multiple of characteristic.
    kgtp_rac: f32,
    /// Rac1 auto-activation rate as a multiple of baseline Rac1 activation rate.
    kgtp_rac_auto: f32,
    /// Chemoattractant factor affecting Rac1 activation rate.
    chemoa: f32,
    /// COA factor affecting Rac1 activation rate.
    coa_half_d: Length,
    /// Baseline Rac1 inactivation rate.
    kdgtp_rac: f32,
    /// RhoA mediated inhibition of Rac1 as a multiple of baseline Rac1 inactivation rate.
    kdgtp_rho_on_rac: f32,
    /// Strain at which Rac1 tension-mediated inhibition is half-strength.
    halfmax_tension_inhib: f32,
    /// Tension-mediated Rac1 inhibition as a multiple of baseline Rac1 inactivation rate.
    tension_inhib: f32,
    /// Baseline RhoA activation rate.
    kgtp_rho: f32,
    /// RhoA auto-activation rate as a multiple of baseline RhoA activation rate.
    kgtp_auto_rho: f32,
    /// Baseline RhoA inactivation rate.
    kdgtp_rho: f32,
    /// Rac1 mediated inhibition of RhoA as a multiple of baseline RhoA inactivation rate.
    kdgtp_rac_on_rho: f32,
    /// CIL factor affecting RhoA activation rate.
    cil: f32,
    rand_scheme: RandomizationType,
    /// Average period between randomization events.
    rand_avg_t: Time,
    /// Standard deviation of period between period of randomization events.
    rand_std_t: Time,
    /// Randomization factors affecting Rac1 activation rate.
    rand: f32,
    /// Fraction of vertices to be selected for increased Rac1 activation due to random events.
    rand_vs: f32,
}

#[derive(Clone)]
pub struct Parameters {
    /// Cell radius.
    pub cell_r: f32,
    /// Resting edge length.
    pub rest_edge_len: f32,
    /// Resting area.
    pub rest_area: f32,
    /// Length criterion to determine if two vertices are "close".
    pub close_criterion: f32,
    /// Viscosity.
    pub vertex_eta: f32,
    /// Rho GTPase at max force.
    pub vertex_rgtp_act_at_max_f: f32,
    /// Stiffness of edge.
    pub stiffness_edge: f32,
    /// RhoA mediated retractive force constant.
    pub max_retractive_f: f32,
    /// Max protrusive force.
    pub max_protrusive_f: f32,
    /// Stiffness of cytoplasm.
    pub stiffness_ctyo: f32,
    /// Rate of Rho GTPase GDI unbinding and subsequent membrane attachment.
    pub k_mem_on: f32,
    /// Rate of Rho GTPase membrane disassociation.
    pub k_mem_off: f32,
    /// Diffusion rate of Rho GTPase on membrane.
    pub diffusion_rgtp: f32,
    /// Initial fraction of Rho GTPase in cytosol (inactive).
    pub init_frac_cyto: f32,
    /// Initial fraction of Rho GTPase active.
    pub init_frac_active: f32,
    /// Initial fraction of Rho GTPase inactive.
    pub init_frac_inactive: f32,
    /// Hill function exponent for Rho GTPase chemistry.
    pub h_exp: f32,
    /// Halfmax Rho GTPase activity.
    pub halfmax_vertex_rgtp_conc: f32,
    /// Total amount of Rac1 in cell.
    pub tot_rac: f32,
    /// Total amount of RhoA in cell.
    pub tot_rho: f32,
    /// Baseline Rac1 activation rate.
    pub kgtp_rac: f32,
    /// Rac1 auto-activation rate as a multiple of baseline Rac1 activation rate.
    pub kgtp_rac_auto: f32,
    /// Chemoattractant factor affecting Rac1 activation rate.
    pub chemoa: f32,
    /// Baseline Rac1 inactivation rate.
    pub kdgtp_rac: f32,
    /// RhoA mediated inhibition of Rac1 as a multiple of baseline Rac1 inactivation rate.
    pub kdgtp_rho_on_rac: f32,
    /// Strain at which Rac1 tension-mediated inhibition is half-strength.
    pub halfmax_tension_inhib: f32,
    /// Tension-mediated Rac1 inhibition as a multiple of baseline Rac1 inactivation rate.
    pub tension_inhib: f32,
    /// Baseline RhoA activation rate.
    pub kgtp_rho: f32,
    /// RhoA auto-activation rate as a multiple of baseline RhoA activation rate.
    pub kgtp_rho_auto: f32,
    /// Baseline RhoA inactivation rate.
    pub kdgtp_rho: f32,
    /// Rac1 mediated inhibition of RhoA as a multiple of baseline RhoA inactivation rate.
    pub kdgtp_rac_on_rho: f32,
    /// CIL factor affecting RhoA activation rate.
    pub cil: f32,
    pub rand_scheme: RandomizationType,
    /// Average period of randomization events.
    pub rand_avg_t: f32,
    /// Standard deviation period of randomization events.
    pub rand_std_t: f32,
    /// Randomization factors affecting Rac1 activation rate.
    pub rand: f32,
    /// Fraction of vertices to be selected for increased Rac1 activation due to random events.
    pub rand_vs: f32,
}

impl Default for UserParams {
    fn default() -> Self {
        let rgtp_d = (Length(0.1).micro().pow(2.0).g() / Time(1.0).g())
            .to_diffusion()
            .unwrap();

        UserParams {
            cell_diam: Length(40.0).micro(),
            close_criterion: Length(0.5).micro(),
            rgtp_act_at_max_f: 0.4,
            stiffness_cortex: Stress(8.0).kilo(),
            lm_h: Length(300.0).nano(),
            halfmax_rgtp_frac: 0.4,
            lm_ss: Stress(10.0).kilo(),
            rho_friction: 0.2,
            stiffness_ctyo: Force(1e-5),
            diffusion_rgtp: rgtp_d,
            init_frac_cyto: 0.8,
            init_frac_active: 0.1,
            h_exp: 3.0,
            tot_rac: 2.5e6,
            tot_rho: 1e6,
            kgtp_rac: 24.0,
            kgtp_rac_auto: 500.0,
            chemoa: 7.5,
            coa_half_d: Length(110.0e-6),
            kdgtp_rac: 4.0,
            kdgtp_rho_on_rac: 4000.0,
            halfmax_tension_inhib: 0.1,
            tension_inhib: 40.0,
            kgtp_rho: 28.0,
            kgtp_auto_rho: 390.0,
            kdgtp_rho: 60.0,
            kdgtp_rac_on_rho: 400.0,
            cil: 60.0,
            rand_scheme: RandomizationType::Off,
            rand_avg_t: Time(40.0 * 60.0),
            rand_std_t: Time(10.0 * 60.0),
            rand: 10.0,
            rand_vs: 0.25,
        }
    }
}

impl UserParams {
    pub fn gen_params(&self, cq: &CharQuants) -> Parameters {
        let cell_r = self.cell_diam.mulf(0.5);
        let rel = self.cell_diam.mulf((PI / (NVERTS as f32)).sin());
        let ra = cell_r.pow(2.0).mulf(PI);
        let close_criterion = self.close_criterion.pow(2.0);
        let max_protrusion_f = self.lm_h.g() * self.lm_ss.g() * rel.g();
        let max_retraction_f = max_protrusion_f.mulf(self.rho_friction);
        let halfmax_rgtp_conc = rel.pow(-1.0).mulf(0.4 / (NVERTS as f32));

        Parameters {
            cell_r: cq.normalize(&cell_r),
            rest_edge_len: cq.normalize(&rel),
            rest_area: cq.normalize(&ra),
            close_criterion: cq.normalize(&close_criterion),
            vertex_eta: cq.normalize(&cq.eta) / (NVERTS as f32),
            vertex_rgtp_act_at_max_f: self.rgtp_act_at_max_f / (NVERTS as f32 / 2.0),
            stiffness_edge: cq.normalize(&(self.stiffness_cortex.g() * cq.l3d.g())),
            max_retractive_f: cq.normalize(&max_retraction_f),
            max_protrusive_f: cq.normalize(&max_protrusion_f),
            stiffness_ctyo: cq.normalize(&self.stiffness_ctyo),
            k_mem_on: cq.normalize(&cq.k_mem_on),
            k_mem_off: cq.normalize(&cq.k_mem_off),
            diffusion_rgtp: cq.normalize(&self.diffusion_rgtp),
            init_frac_cyto: self.init_frac_active,
            init_frac_active: self.init_frac_cyto,
            init_frac_inactive: {
                let ifi = 1.0 - self.init_frac_active - self.init_frac_cyto;
                if !(ifi < 0.0) {
                    ifi
                } else {
                    panic!("Initial fraction of inactive Rho GTPase is negative.")
                }
            },
            h_exp: self.h_exp,
            halfmax_vertex_rgtp_conc: cq.normalize(&halfmax_rgtp_conc) / (NVERTS as f32 / 2.0),
            tot_rac: self.tot_rac,
            tot_rho: self.tot_rho,
            kgtp_rac: cq.normalize(&cq.kgtp.mulf(self.kgtp_rac)),
            kgtp_rac_auto: cq.normalize(&cq.kgtp.mulf(self.kgtp_rac_auto)),
            chemoa: self.chemoa,
            kdgtp_rac: cq.normalize(&cq.kdgtp.mulf(self.kdgtp_rac)),
            kdgtp_rho_on_rac: cq.normalize(&cq.kdgtp.mulf(self.kdgtp_rho_on_rac)),
            halfmax_tension_inhib: self.halfmax_tension_inhib,
            tension_inhib: self.tension_inhib,
            kgtp_rho: cq.normalize(&cq.kgtp.mulf(self.kgtp_rho)),
            kgtp_rho_auto: cq.normalize(&cq.kgtp.mulf(self.kgtp_auto_rho)),
            kdgtp_rho: cq.normalize(&cq.kdgtp.mulf(self.kdgtp_rho)),
            kdgtp_rac_on_rho: cq.normalize(&cq.kdgtp.mulf(self.kdgtp_rac_on_rho)),
            cil: self.cil,
            rand_scheme: self.rand_scheme.clone(),
            rand_avg_t: cq.normalize(&self.rand_avg_t),
            rand_std_t: cq.normalize(&self.rand_std_t),
            rand: self.rand,
            rand_vs: self.rand_vs,
        }
    }
}
