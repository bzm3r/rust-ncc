// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod quantity;

use crate::cell::calc_init_cell_area;
use crate::cell::chemistry::RgtpDistribution;
use crate::interactions::CilMat;
use crate::parameters::quantity::{
    Diffusion, Force, Length, Quantity, Stress, Time, Tinv, Viscosity,
};
use crate::NVERTS;
use std::f32::consts::PI;

/// Characteristic quantities used for normalization.
pub struct BasicQuants {
    pub eta: Viscosity,
    pub f: Force,
    pub l: Length,
    pub t: Time,
    pub l3d: Length,
    pub k_mem_off: Tinv,
    pub k_mem_on: Tinv,
    pub kgtp: Tinv,
    pub kdgtp: Tinv,
    pub frac_rgtp: f32,
}

impl BasicQuants {
    pub fn normalize<T: Quantity>(&self, q: &T) -> f32 {
        let q = q.g();
        let u = q.units();
        (q * self.f.pow(-1.0 * u.f) * self.l.pow(-1.0 * u.l) * self.t.pow(-1.0 * u.t)).value()
    }

    pub fn time(&self) -> f32 {
        self.t.0
    }
}

pub struct RawWorldParameters {
    pub vertex_eta: Viscosity,
    pub close_criterion: Length,
    pub cil: CilMat,
    pub adh_const: Force,
}

#[derive(Clone)]
pub struct WorldParameters {
    pub vertex_eta: f32,
    pub close_criterion: f32,
    pub cil: CilMat,
    pub adh_const: f32,
}

impl RawWorldParameters {
    pub fn normalize(&self, bq: &BasicQuants) -> WorldParameters {
        WorldParameters {
            vertex_eta: bq.normalize(&self.vertex_eta),
            close_criterion: bq.normalize(&self.close_criterion),
            cil: self.cil.clone(),
            adh_const: bq.normalize(&self.adh_const),
        }
    }
}

pub struct RawParameters {
    /// Cell diameter.
    pub cell_diam: Length,
    /// Fraction of max force achieved at `rgtp_act_at_max_f`.
    pub halfmax_rgtp_max_f_frac: f32,
    /// Stiffness of the membrane-cortex complex.
    pub stiffness_cortex: Stress,
    /// Typical lamellipod height.
    pub lm_h: Length,
    /// Halfmax Rho GTPase activity.
    pub halfmax_rgtp_frac: f32,
    /// Lamellipod stall stress.
    pub lm_ss: Stress,
    /// Friction force opposing RhoA pulling.
    pub rho_friction: f32,
    /// Stiffness of cytoplasm.
    pub stiffness_ctyo: Force,
    /// Diffusion rate of Rho GTPase on membrane as a multiple of characteristic.
    pub diffusion_rgtp: Diffusion,
    /// Initial distribution of Rac1.
    pub init_rac: RgtpDistribution,
    /// Initial distribution of RhoA.
    pub init_rho: RgtpDistribution,
    /// Total amount of Rac1 in cell.
    pub tot_rac: f32,
    /// Total amount of RhoA in cell.
    pub tot_rho: f32,
    /// Baseline Rac1 activation rate as a multiple of characteristic.
    pub kgtp_rac: f32,
    /// Rac1 auto-activation rate as a multiple of baseline Rac1 activation rate.
    pub kgtp_rac_auto: f32,
    /// Chemoattractant factor affecting Rac1 activation rate.
    pub chemoa: f32,
    /// COA factor affecting Rac1 activation rate.
    pub coa_half_d: Length,
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
    pub kgtp_auto_rho: f32,
    /// Baseline RhoA inactivation rate.
    pub kdgtp_rho: f32,
    /// Rac1 mediated inhibition of RhoA as a multiple of baseline RhoA inactivation rate.
    pub kdgtp_rac_on_rho: f32,
    /// Enable randomization?
    pub randomization: bool,
    /// Average period between randomization events.
    pub rand_avg_t: Time,
    /// Magnitude of randomly applied factor affecting Rac1 activation rate.
    pub rand_mag: f32,
    /// Fraction of vertices to be selected for increased Rac1 activation due to random events.
    pub rand_vs: f32,
    pub rand_std_t: Time,
}

#[derive(Clone)]
pub struct Parameters {
    /// Cell radius.
    pub cell_r: f32,
    /// Resting edge length.
    pub rest_edge_len: f32,
    /// Resting area.
    pub rest_area: f32,
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
    /// Initial distribution of Rac1.
    pub init_rac: RgtpDistribution,
    /// Initial distribution of RhoA.
    pub init_rho: RgtpDistribution,
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
    /// Enable randomization?
    pub randomization: bool,
    /// Average time between random events, in timesteps.
    pub rand_avg_t: f32,
    /// Standard deviation of time between random events, in timesteps.
    pub rand_std_t: f32,
    /// Magnitude of factor randomly applied to Rac1 activation rate.
    pub rand_mag: f32,
    /// Number of vertices to be selected for random Rac1 activity boost.
    pub num_rand_vs: usize,
    pub total_rgtp: f32,
}

impl RawParameters {
    pub fn gen_parameters(&self, bq: &BasicQuants) -> Parameters {
        let cell_r = self.cell_diam.mulf(0.5);
        let rel = self.cell_diam.mulf((PI / (NVERTS as f32)).sin());
        let ra = Length(1.0)
            .pow(2.0)
            .mulf(calc_init_cell_area(cell_r.value(), NVERTS));
        let const_protrusive =
            (self.lm_h.g() * self.lm_ss.g() * rel.g()).mulf(self.halfmax_rgtp_max_f_frac);
        let const_retractive = const_protrusive.mulf(self.rho_friction);
        let halfmax_vertex_rgtp_act = (self.halfmax_rgtp_frac / bq.frac_rgtp) / NVERTS as f32;
        let halfmax_vertex_rgtp_conc = rel.pow(-1.0).mulf(halfmax_vertex_rgtp_act);
        let stiffness_edge = self.stiffness_cortex.g() * bq.l3d.g();
        let stiffness_cyto = self.stiffness_ctyo.g().mulf(1.0 / NVERTS as f32);

        Parameters {
            cell_r: bq.normalize(&cell_r),
            rest_edge_len: bq.normalize(&rel),
            rest_area: bq.normalize(&ra),
            stiffness_edge: bq.normalize(&stiffness_edge),
            const_protrusive: bq.normalize(&const_protrusive),
            const_retractive: bq.normalize(&const_retractive),
            stiffness_ctyo: bq.normalize(&stiffness_cyto),
            k_mem_on_vertex: bq.normalize(&bq.k_mem_on) / NVERTS as f32,
            k_mem_off: bq.normalize(&bq.k_mem_off),
            diffusion_rgtp: bq.normalize(&self.diffusion_rgtp),
            init_rac: self.init_rac,
            init_rho: self.init_rho,
            halfmax_vertex_rgtp_act,
            halfmax_vertex_rgtp_conc: bq.normalize(&halfmax_vertex_rgtp_conc),
            tot_rac: self.tot_rac,
            tot_rho: self.tot_rho,
            kgtp_rac: bq.normalize(&bq.kgtp.mulf(self.kgtp_rac)),
            kgtp_rac_auto: bq.normalize(&bq.kgtp.mulf(self.kgtp_rac_auto)),
            chemoa: self.chemoa,
            kdgtp_rac: bq.normalize(&bq.kdgtp.mulf(self.kdgtp_rac)),
            kdgtp_rho_on_rac: bq.normalize(&bq.kdgtp.mulf(self.kdgtp_rho_on_rac)),
            halfmax_tension_inhib: self.halfmax_tension_inhib,
            tension_inhib: self.tension_inhib,
            kgtp_rho: bq.normalize(&bq.kgtp.mulf(self.kgtp_rho)),
            kgtp_rho_auto: bq.normalize(&bq.kgtp.mulf(self.kgtp_auto_rho)),
            kdgtp_rho: bq.normalize(&bq.kdgtp.mulf(self.kdgtp_rho)),
            kdgtp_rac_on_rho: bq.normalize(&bq.kdgtp.mulf(self.kdgtp_rac_on_rho)),
            randomization: self.randomization,
            rand_avg_t: bq.normalize(&self.rand_avg_t).ceil(),
            rand_std_t: bq.normalize(&self.rand_std_t).ceil(),
            rand_mag: self.rand_mag,
            num_rand_vs: (self.rand_vs * NVERTS as f32) as usize,
            total_rgtp: 1.0 / bq.frac_rgtp,
        }
    }
}
