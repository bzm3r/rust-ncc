// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::consts::NVERTS;
use crate::quantity::{Diffusion, Force, Length, Quantity, Stress, Time, Tinv, Viscosity};
use crate::random::RandomizationType;
use override_derive::Overrides;
use std::f32::consts::PI;
use crate::math::calc_init_cell_area;

/// World-wide parameters.
pub struct WorldParameters {
    eta: Viscosity,
    f: Force,
    l: Length,
    t: Time,
    l3d: Length,
    k_mem_off: Tinv,
    k_mem_on: Tinv,
    kgtp: Tinv,
    kdgtp: Tinv,
}

impl WorldParameters {
    pub fn normalize<T: Quantity>(&self, q: &T) -> f32 {
        let q = q.g();
        let u = q.units();
        (q * self.f.pow(-1.0 * u.f) * self.l.pow(-1.0 * u.l) * self.t.pow(-1.0 * u.t)).value()
    }

    pub fn time(&self) -> f32 {
        self.t.0
    }
}

impl Default for WorldParameters {
    fn default() -> Self {
        // Stress on lamellipod is on order of 1kPa, height of lamellipod on order of 100 nm, length of edge on order of 10 um
        let f = (Stress(1.0).kilo().g() * Length(100.0).nano().g() * Length(10.0).micro().g())
            .to_force()
            .unwrap();
        let eta = Viscosity(0.29);

        WorldParameters {
            eta,
            l: Length(1.0).micro(),
            t: Time(2.0),
            f,
            l3d: Length(10e-6),
            k_mem_off: Tinv(0.15),
            k_mem_on: Tinv(0.02),
            kgtp: Tinv(1e-4),
            kdgtp: Tinv(1e-4),
            //coa_dist: Length(110.0).micro(),
        }
    }
}

#[derive(Overrides)]
pub struct InputParameters {
    /// Cell diameter.
    pub cell_diam: Length,
    /// Length criterion to determine if two vertices are "close".
    close_criterion: Length,
    /// Fraction of max force achieved at `rgtp_act_at_max_f`.
    halfmax_rgtp_max_f_frac: f32,
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
    init_frac_inactive: f32,
    /// Initial fraction of Rho GTPase active.
    init_frac_active: f32,
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
    /// Stiffness of edge.
    pub stiffness_edge: f32,
    /// Rac1 mediated protrusive force constant.
    pub const_protrusive: f32,
    /// RhoA mediated protrusive force constant.
    pub const_retractive: f32,
    /// Stiffness of cytoplasm.
    pub stiffness_ctyo: f32,
    /// Rate of Rho GTPase GDI unbinding and subsequent membrane attachment.
    pub k_mem_on_vertex: f32,
    /// Rate of Rho GTPase membrane disassociation.
    pub k_mem_off: f32,
    /// Diffusion rate of Rho GTPase on membrane.
    pub diffusion_rgtp: f32,
    /// Initial fraction of Rho GTPase inactive.
    pub init_frac_inactive: f32,
    /// Initial fraction of Rho GTPase active.
    pub init_frac_active: f32,
    /// Halfmax Rho GTPase activity per vertex.
    pub halfmax_vertex_rgtp_act: f32,
    /// Halfmax Rho GTPase activity per vertex as concentration.
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

impl Default for InputParameters {
    fn default() -> Self {
        let rgtp_d = (Length(0.1_f32.sqrt()).micro().pow(2.0).g() / Time(1.0).g())
            .to_diffusion()
            .unwrap();

        InputParameters {
            cell_diam: Length(40.0).micro(),
            close_criterion: Length(0.5).micro(),
            stiffness_cortex: Stress(8.0).kilo(),
            lm_h: Length(200.0).nano(),
            halfmax_rgtp_max_f_frac: 0.3,
            halfmax_rgtp_frac: 0.4,
            lm_ss: Stress(10.0).kilo(),
            rho_friction: 0.2,
            stiffness_ctyo: Force(1e-5),
            diffusion_rgtp: rgtp_d,
            init_frac_inactive: 0.1,
            init_frac_active: 0.1,
            tot_rac: 2.5e6,
            tot_rho: 1e6,
            kgtp_rac: 24.0,
            kgtp_rac_auto: 500.0,
            chemoa: 7.5,
            coa_half_d: Length(110.0e-6),
            kdgtp_rac: 8.0,
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

impl InputParameters {
    pub fn gen_parameters(&self, cq: &WorldParameters) -> Parameters {
        let cell_r = self.cell_diam.mulf(0.5);
        let rel = self.cell_diam.mulf((PI / (NVERTS as f32)).sin());
        let ra = Length(1.0).pow(2.0).mulf(calc_init_cell_area(cell_r.value(), NVERTS));
        let close_criterion = self.close_criterion.pow(2.0);
        let const_protrusive =
            (self.lm_h.g() * self.lm_ss.g() * rel.g()).mulf(self.halfmax_rgtp_max_f_frac);
        let const_retractive = const_protrusive.mulf(self.rho_friction);
        let halfmax_vertex_rgtp_act = self.halfmax_rgtp_frac / NVERTS as f32;
        let halfmax_vertex_rgtp_conc = rel.pow(-1.0).mulf(halfmax_vertex_rgtp_act);
        let stiffness_edge = self.stiffness_cortex.g() * cq.l3d.g();
        let stiffness_cyto = self.stiffness_ctyo.g().mulf(1.0 / NVERTS as f32);
        let (init_frac_active, init_frac_inactive) = if 1.0 - self.init_frac_active - self.init_frac_inactive < 0.0 {
            panic!("Cytosolic fraction is negative. init_frac_active: {}, init_frac_inactive: {}", self.init_frac_active, self.init_frac_inactive);
        } else {
            (self.init_frac_active, self.init_frac_inactive)
        };
        Parameters {
            cell_r: cq.normalize(&cell_r),
            rest_edge_len: cq.normalize(&rel),
            rest_area: cq.normalize(&ra),
            close_criterion: cq.normalize(&close_criterion),
            vertex_eta: cq.normalize(&cq.eta) / (NVERTS as f32),
            stiffness_edge: cq.normalize(&stiffness_edge),
            const_protrusive: cq.normalize(&const_protrusive),
            const_retractive: cq.normalize(&const_retractive),
            stiffness_ctyo: cq.normalize(&stiffness_cyto),
            k_mem_on_vertex: cq.normalize(&cq.k_mem_on) / NVERTS as f32,
            k_mem_off: cq.normalize(&cq.k_mem_off),
            diffusion_rgtp: cq.normalize(&self.diffusion_rgtp),
            init_frac_inactive,
            init_frac_active,
            halfmax_vertex_rgtp_act,
            halfmax_vertex_rgtp_conc: cq.normalize(&halfmax_vertex_rgtp_conc),
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
