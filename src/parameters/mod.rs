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
use crate::math::geometry::BBox;
use crate::math::v2d::V2D;
use crate::parameters::quantity::{
    Diffusion, Force, Length, Quantity, Stress, Time, Tinv, Viscosity,
};
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::f32::consts::PI;

/// Characteristic quantities used for normalization.
#[derive(Clone, Copy, Deserialize, Serialize, Default, Debug)]
pub struct CharQuantities {
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

impl CharQuantities {
    /// Given a quantity `q`, normalize its units using the primary units `f` (Force),
    /// `l` (`Length`) and `t` (`Time`) provided in `CharQuants`.
    pub fn normalize<T: Quantity>(&self, q: &T) -> f32 {
        let q = q.g();
        let u = q.units();
        (q * self.f.pow(-1.0 * u.f)
            * self.l.pow(-1.0 * u.l)
            * self.t.pow(-1.0 * u.t))
        .number()
    }

    pub fn time(&self) -> f32 {
        self.t.0
    }
}

#[derive(Clone, Debug)]
pub struct RawCloseBounds {
    pub zero_at: Length,
    pub one_at: Length,
}

impl RawCloseBounds {
    pub fn new(zero_at: Length, one_at: Length) -> RawCloseBounds {
        RawCloseBounds { zero_at, one_at }
    }
}

#[derive(Clone, Debug)]
pub struct RawPhysicalContactParams {
    pub range: RawCloseBounds,
    pub adh_mag: Option<Force>,
    pub cal_mag: Option<f32>,
    pub cil_mag: f32,
}

impl RawPhysicalContactParams {
    pub fn refine(
        &self,
        cq: &CharQuantities,
    ) -> PhysicalContactParams {
        PhysicalContactParams {
            range: CloseBounds::new(
                cq.normalize(&self.range.zero_at),
                cq.normalize(&self.range.one_at),
            ),
            adh_mag: self
                .adh_mag
                .map(|adh_mag| cq.normalize(&adh_mag)),
            cal_mag: self.cal_mag,
            cil_mag: self.cil_mag,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RawCoaParams {
    /// Factor controlling to what extent line-of-sight blockage should be penalized.
    pub los_penalty: f32,
    pub range: Length,
    pub mag: f32,
}

impl RawCoaParams {
    pub fn refine(&self, bq: &CharQuantities) -> CoaParams {
        let range = bq.normalize(&self.range);
        CoaParams {
            los_penalty: self.los_penalty,
            range,
            mag: self.mag,
            // self.mag * exp(distrib_exp * x), where x is distance
            // between points.
            distrib_exp: 0.5f32.ln() / (0.5 * range),
        }
    }
}

#[derive(Clone, Debug)]
pub struct RawChemAttrParams {
    center: [Length; 2],
    center_mag: f32,
    drop_per_char_l: f32,
    char_l: Length,
}

impl RawChemAttrParams {
    pub fn refine(&self, bq: &CharQuantities) -> ChemAttrParams {
        ChemAttrParams {
            center: V2D {
                x: bq.normalize(&self.center[0]),
                y: bq.normalize(&self.center[1]),
            },
            center_mag: self.center_mag,
            slope: self.drop_per_char_l / bq.normalize(&self.char_l),
        }
    }
}

#[derive(Clone, Debug)]
pub struct RawBdryParams {
    shape: Vec<[Length; 2]>,
    skip_bb_check: bool,
    mag: f32,
}

impl RawBdryParams {
    pub fn refine(&self, bq: &CharQuantities) -> BdryParams {
        let shape = self
            .shape
            .iter()
            .map(|p| V2D {
                x: bq.normalize(&p[0]),
                y: bq.normalize(&p[1]),
            })
            .collect::<Vec<V2D>>();
        let bbox = BBox::from_points(&shape);
        BdryParams {
            shape,
            bbox,
            skip_bb_check: self.skip_bb_check,
            mag: self.mag,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RawInteractionParams {
    pub coa: Option<RawCoaParams>,
    pub chem_attr: Option<RawChemAttrParams>,
    pub bdry: Option<RawBdryParams>,
    pub phys_contact: RawPhysicalContactParams,
}

impl RawInteractionParams {
    pub fn refine(&self, bq: &CharQuantities) -> InteractionParams {
        InteractionParams {
            coa: self.coa.as_ref().map(|coa| coa.refine(bq)),
            chem_attr: self
                .chem_attr
                .as_ref()
                .map(|chem_attr| chem_attr.refine(bq)),
            bdry: self.bdry.as_ref().map(|bdry| bdry.refine(bq)),
            phys_contact: self.phys_contact.refine(bq),
        }
    }
}

#[derive(Clone, Debug)]
pub struct RawWorldParameters {
    pub vertex_eta: Viscosity,
    pub interactions: RawInteractionParams,
}

#[derive(
    Clone, Copy, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct CloseBounds {
    pub zero_at: f32,
    pub one_at: f32,
}

impl CloseBounds {
    pub fn new(zero_until: f32, one_at: f32) -> CloseBounds {
        CloseBounds {
            zero_at: zero_until,
            one_at,
        }
    }
}

#[derive(
    Clone, Copy, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct PhysicalContactParams {
    /// Maximum distance between two points, for them to be considered
    /// in contact. This is usually set to 0.5 micrometers.
    pub range: CloseBounds,
    /// Optional adhesion magnitude. If it is `None`, no adhesion
    /// will be calculated.
    pub adh_mag: Option<f32>,
    /// Optional CAL magnitude. If it is `None`, simulation will
    /// always execute CIL upon contact.
    pub cal_mag: Option<f32>,
    /// Magnitude of CIL that acts on Rho GTPase activation/
    /// inactivation rates.
    pub cil_mag: f32,
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub struct CoaParams {
    //TODO: Expand upon LOS system.
    /// Factor controlling to what extent line-of-sight blockage
    /// should be penalized. See SI for further information.
    pub los_penalty: f32,
    //TODO: Expand upon the exponential curve used to mdoel COA
    // (see SI). Confirm whether it is half-max, or full-max range.
    /// Factor controlling shape of the exponential modelling COA
    /// interaction. It captures the (half?)-maximum distance between
    /// two points still able to undergo COA.
    pub range: f32,
    /// Magnitude of COA that acts on Rac1 activation rates.
    pub mag: f32,
    //TODO: look up exactly what is being done for this (see where
    // parameter is being generated for hint).
    /// Factor controlling the shape of the exponential modelling
    /// COA interaction (a function shaping parameter). It determines
    /// the distance at which two points would sense COA at half-max
    /// magnitude.
    pub distrib_exp: f32,
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub struct ChemAttrParams {
    /// Location of the chemoattractant center.
    pub center: V2D,
    /// Magnitude of chemoattractant a cell would sense if it were
    /// right on top of the chemoattractant source.
    pub center_mag: f32,
    /// Assuming shallow chemoattractant gradient, which can be
    /// modelled using a linear function with slope `slope`.
    pub slope: f32,
}

#[derive(Clone, Deserialize, Serialize, PartialEq, Debug)]
pub struct BdryParams {
    /// Shape of the boundary.
    pub shape: Vec<V2D>,
    /// Bounding box of the boundary.
    pub bbox: BBox,
    /// Should boundary bounding box be checked to see if cell is
    /// within the boundary?
    pub skip_bb_check: bool,
    /// Magnitude of CIL-type interaction.
    pub mag: f32,
}

#[derive(Clone, Deserialize, Serialize, PartialEq, Default, Debug)]
pub struct InteractionParams {
    pub phys_contact: PhysicalContactParams,
    pub coa: Option<CoaParams>,
    pub chem_attr: Option<ChemAttrParams>,
    pub bdry: Option<BdryParams>,
}

#[derive(Clone, Deserialize, Serialize, PartialEq, Default, Debug)]
pub struct WorldParameters {
    /// Viscosity value used to calculate change in position of a
    /// vertex due to calculated forces on it.
    pub vertex_eta: f32,
    pub interactions: InteractionParams,
}

impl RawWorldParameters {
    pub fn refine(&self, bq: &CharQuantities) -> WorldParameters {
        WorldParameters {
            vertex_eta: bq.normalize(&self.vertex_eta),
            interactions: self.interactions.refine(bq),
        }
    }
}

/// The "raw", unprocessed, parameters that are supplied by the user.
pub struct RawParameters {
    /// Cell diameter.
    pub cell_diam: Length,
    /// Fraction of max force achieved at `rgtp_act_at_max_f`.
    pub halfmax_rgtp_max_f_frac: f32,
    /// Stiffness of the membrane-cortex complex.
    pub stiffness_cortex: Stress,
    /// Typical lamellipod height: typical height of lamellipod (on the order of 100 nm).
    pub lm_h: Length,
    /// Halfmax Rho GTPase activity.
    pub halfmax_rgtp_frac: f32,
    /// Lamellipod stall stress: how much stress can lamellipod exert at most.
    pub lm_ss: Stress,
    /// Friction force opposing RhoA pulling.
    pub rho_friction: f32,
    /// Stiffness of cytoplasm.
    pub stiffness_ctyo: Force,
    /// Diffusion rate of Rho GTPase on membrane.
    pub diffusion_rgtp: Diffusion,
    /// Initial distribution of Rac1.
    pub init_rac: RgtpDistribution,
    /// Initial distribution of RhoA.
    pub init_rho: RgtpDistribution,
    /// Total amount of Rac1 in cell.
    pub tot_rac: f32,
    /// Total amount of RhoA in cell.
    pub tot_rho: f32,
    /// Baseline Rac1 activation rate as a multiple of characteristic activation rate, which we
    /// currently take to be 1e-4/second.
    ///
    /// The user only needs to think about providing the baseline Rac1 activity rate in terms of
    /// the characteristic activity rate.
    pub kgtp_rac: f32,
    /// Rac1 auto-activation rate as a multiple of the characteristic activity rate.
    pub kgtp_rac_auto: f32,
    /// Baseline Rac1 inactivation rate as a multiple of the characteristic inactivity rate.
    pub kdgtp_rac: f32,
    /// RhoA mediated inhibition of Rac1 as a multiple of baseline Rac1 inactivation rate.
    pub kdgtp_rho_on_rac: f32,
    /// Strain at which Rac1 tension-mediated inhibition is half-strength.
    pub halfmax_tension_inhib: f32,
    /// Maximum tension-mediated Rac1 inhibition as a multiple of baseline Rac1 inactivation rate.
    pub tension_inhib: f32,
    /// Baseline RhoA activation rate as a multiple of characteristic activation rate.
    pub kgtp_rho: f32,
    /// RhoA auto-activation rate as a multiple of characteristic activation rate.
    pub kgtp_auto_rho: f32,
    /// Baseline RhoA inactivation rate as a multiple of characteristic inactivation rate.
    pub kdgtp_rho: f32,
    /// Rac1 mediated inhibition of RhoA as a multiple of characteristic inactivation rate.
    pub kdgtp_rac_on_rho: f32,
    /// Enable randomization of bursts in Rac1 activity?
    pub randomization: bool,
    /// Average period between randomization events.
    pub rand_avg_t: Time,
    /// Standard deviation of period between randomization events.
    pub rand_std_t: Time,
    /// Magnitude of randomly applied factor affecting Rac1 activation rate: how big a burst?
    pub rand_mag: f32,
    /// Fraction of vertices to be selected for increased Rac1 activation due to random events.
    pub rand_vs: f32,
}

#[derive(Clone, Deserialize, Serialize, Default, Debug, PartialEq)]
pub struct Parameters {
    /// Resting cell radius.
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
    /// Enable randomization of bursts in Rac1 activity?
    pub randomization: bool,
    /// Average time between random events, in timesteps.
    pub rand_avg_t: f32,
    /// Standard deviation of time between random events, in timesteps.
    pub rand_std_t: f32,
    /// Magnitude of factor randomly applied to Rac1 activation rate.
    pub rand_mag: f32,
    /// Number of vertices to be selected for random Rac1 activity boost.
    pub num_rand_vs: u32,
    /// Total Rho GTPase "fraction" compared to "baseline" fraction of Rho GTPase activity.
    /// 0.05 = "halfmax level of Rho GTPase"
    /// 1.0 = " total amount"
    ///
    /// At a vertex, assume that we are only at 5% of the halfmax level:
    /// 0.05 * 0.05 = 0.0025
    ///
    /// 1.0 = "halfmax level of Rho GTPase"
    /// 20.0 = "total amount"
    /// At a vertex, assume that we are only at 5% of the halfmax:
    /// 0.05 * 1.0 = 0.05
    //TODO: should this be removed?
    pub total_rgtp: f32,
}

impl RawParameters {
    pub fn gen_parameters(&self, bq: &CharQuantities) -> Parameters {
        let cell_r = self.cell_diam.mul_number(0.5);
        let rel =
            self.cell_diam.mul_number((PI / (NVERTS as f32)).sin());
        let ra = Length(1.0)
            .pow(2.0)
            .mul_number(calc_init_cell_area(cell_r.number()));
        let const_protrusive =
            (self.lm_h.g() * self.lm_ss.g() * rel.g())
                .mul_number(self.halfmax_rgtp_max_f_frac);
        let const_retractive =
            const_protrusive.mul_number(self.rho_friction);
        let halfmax_vertex_rgtp_act =
            (self.halfmax_rgtp_frac / bq.frac_rgtp) / NVERTS as f32;
        let halfmax_vertex_rgtp_conc =
            rel.pow(-1.0).mul_number(halfmax_vertex_rgtp_act);
        let stiffness_edge = self.stiffness_cortex.g() * bq.l3d.g();
        let stiffness_cyto =
            self.stiffness_ctyo.g().mul_number(1.0 / NVERTS as f32);

        Parameters {
            cell_r: bq.normalize(&cell_r),
            rest_edge_len: bq.normalize(&rel),
            rest_area: bq.normalize(&ra),
            stiffness_edge: bq.normalize(&stiffness_edge),
            const_protrusive: bq.normalize(&const_protrusive),
            const_retractive: bq.normalize(&const_retractive),
            stiffness_ctyo: bq.normalize(&stiffness_cyto),
            k_mem_on_vertex: bq.normalize(&bq.k_mem_on)
                / NVERTS as f32,
            k_mem_off: bq.normalize(&bq.k_mem_off),
            diffusion_rgtp: bq.normalize(&self.diffusion_rgtp),
            init_rac: self.init_rac,
            init_rho: self.init_rho,
            halfmax_vertex_rgtp_act,
            halfmax_vertex_rgtp_conc: bq
                .normalize(&halfmax_vertex_rgtp_conc),
            tot_rac: self.tot_rac,
            tot_rho: self.tot_rho,
            kgtp_rac: bq
                .normalize(&bq.kgtp.mul_number(self.kgtp_rac)),
            kgtp_rac_auto: bq
                .normalize(&bq.kgtp.mul_number(self.kgtp_rac_auto)),
            kdgtp_rac: bq
                .normalize(&bq.kdgtp.mul_number(self.kdgtp_rac)),
            kdgtp_rho_on_rac: bq.normalize(
                &bq.kdgtp.mul_number(self.kdgtp_rho_on_rac),
            ),
            halfmax_tension_inhib: self.halfmax_tension_inhib,
            tension_inhib: self.tension_inhib,
            kgtp_rho: bq
                .normalize(&bq.kgtp.mul_number(self.kgtp_rho)),
            kgtp_rho_auto: bq
                .normalize(&bq.kgtp.mul_number(self.kgtp_auto_rho)),
            kdgtp_rho: bq
                .normalize(&bq.kdgtp.mul_number(self.kdgtp_rho)),
            kdgtp_rac_on_rho: bq.normalize(
                &bq.kdgtp.mul_number(self.kdgtp_rac_on_rho),
            ),
            randomization: self.randomization,
            rand_avg_t: bq.normalize(&self.rand_avg_t).ceil(),
            rand_std_t: bq.normalize(&self.rand_std_t).ceil(),
            rand_mag: self.rand_mag,
            num_rand_vs: (self.rand_vs * NVERTS as f32) as u32,
            total_rgtp: 1.0 / bq.frac_rgtp,
        }
    }
}
