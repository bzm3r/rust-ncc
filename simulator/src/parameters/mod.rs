// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod quantity;
use crate::cell::chemistry::RgtpDistribution;
use crate::math::geometry::{calc_poly_area, BBox};
use crate::math::v2d::V2d;
use crate::parameters::quantity::{
    Diffusion, Force, Length, Quantity, Stress, Time, Tinv, Viscosity,
};
use crate::NVERTS;
use modify_derive::Modify;
use rand_distr::num_traits::Pow;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Characteristic quantities used for normalization.
#[derive(
    Clone,
    Copy,
    Deserialize,
    Serialize,
    Default,
    Debug,
    PartialEq,
    Modify,
)]
pub struct CharQuantities {
    pub eta: Viscosity,
    pub f: Force,
    pub l: Length,
    pub t: Time,
    pub l3d: Length,
    pub kgtp: Tinv,
}

impl CharQuantities {
    /// Given a quantity `q`, normalize its units using the primary units `f` (Force),
    /// `l` (`Length`) and `t` (`Time`) provided in `CharQuants`.
    pub fn normalize<T: Quantity>(&self, q: &T) -> f64 {
        let q = q.g();
        let u = q.units();
        (q * self.f.pow(-1.0 * u.f)
            * self.l.pow(-1.0 * u.l)
            * self.t.pow(-1.0 * u.t))
        .number()
    }

    pub fn time(&self) -> f64 {
        self.t.0
    }
}

#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    Deserialize,
    Serialize,
    PartialEq,
    Modify,
)]
pub struct RawCloseBounds {
    pub zero_at: Length,
    pub one_at: Length,
}

impl RawCloseBounds {
    pub fn new(zero_at: Length, one_at: Length) -> RawCloseBounds {
        RawCloseBounds { zero_at, one_at }
    }
}

#[derive(
    Copy,
    Clone,
    Debug,
    Default,
    Deserialize,
    Serialize,
    PartialEq,
    Modify,
)]
pub struct RawPhysicalContactParams {
    pub range: RawCloseBounds,
    pub adh_mag: Option<Force>,
    pub cal_mag: Option<f64>,
    pub cil_mag: f64,
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

#[derive(
    Deserialize,
    Serialize,
    Clone,
    Copy,
    PartialEq,
    Default,
    Debug,
    Modify,
)]
pub struct RawCoaParams {
    /// Factor controlling to what extent line-of-sight blockage should be
    /// penalized.
    pub los_penalty: f64,
    /// Distance from point of emission at which COA signal reaches half
    /// its maximum value.
    pub halfmax_dist: Length,
    /// Magnitude of COA. It will be divided by `NVERTS` so that it scales based
    /// on the number of vertices.
    pub mag: f64,
    /// If two vertices are within this distance, then COA cannot occur between them.
    pub too_close_dist: Length,
}

impl RawCoaParams {
    pub fn refine(&self, bq: &CharQuantities) -> CoaParams {
        let halfmax_dist = bq.normalize(&self.halfmax_dist);
        CoaParams {
            los_penalty: self.los_penalty,
            halfmax_dist,
            vertex_mag: self.mag / NVERTS as f64,
            // self.mag * exp(distrib_exp * x), where x is distance
            // between points.
            distrib_exp: 0.5f64.ln() / halfmax_dist,
            too_close_dist_sq: bq
                .normalize(&self.too_close_dist)
                .pow(2),
        }
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct RawChemAttrParams {
    center: [Length; 2],
    center_mag: f64,
    drop_per_char_l: f64,
    char_l: Length,
}

impl RawChemAttrParams {
    pub fn refine(&self, bq: &CharQuantities) -> ChemAttrParams {
        ChemAttrParams {
            center: V2d {
                x: bq.normalize(&self.center[0]),
                y: bq.normalize(&self.center[1]),
            },
            center_mag: self.center_mag,
            slope: self.drop_per_char_l / bq.normalize(&self.char_l),
        }
    }
}

#[derive(
    Deserialize, Serialize, Clone, Copy, PartialEq, Default, Debug,
)]
pub struct RawBdryParams {
    shape: [[Length; 2]; 4],
    skip_bb_check: bool,
    mag: f64,
}

impl RawBdryParams {
    pub fn refine(&self, bq: &CharQuantities) -> BdryParams {
        let shape = self
            .shape
            .iter()
            .map(|p| V2d {
                x: bq.normalize(&p[0]),
                y: bq.normalize(&p[1]),
            })
            .collect::<Vec<V2d>>();
        let bbox = BBox::from_points(&shape);
        BdryParams {
            shape,
            bbox,
            skip_bb_check: self.skip_bb_check,
            mag: self.mag,
        }
    }
}

#[derive(
    Deserialize,
    Serialize,
    Clone,
    Copy,
    PartialEq,
    Default,
    Debug,
    Modify,
)]
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

#[derive(
    Deserialize,
    Serialize,
    Copy,
    Clone,
    PartialEq,
    Default,
    Debug,
    Modify,
)]
pub struct RawWorldParameters {
    pub vertex_eta: Viscosity,
    pub interactions: RawInteractionParams,
}

#[derive(
    Clone, Copy, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct CloseBounds {
    pub zero_at_sq: f64,
    pub zero_at: f64,
    pub one_at: f64,
}

impl CloseBounds {
    pub fn new(zero_at: f64, one_at: f64) -> CloseBounds {
        CloseBounds {
            zero_at_sq: zero_at.pow(2),
            zero_at,
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
    pub adh_mag: Option<f64>,
    /// Optional CAL magnitude. If it is `None`, simulation will
    /// always execute CIL upon contact.
    pub cal_mag: Option<f64>,
    /// Magnitude of CIL that acts on Rho GTPase activation/
    /// inactivation rates.
    pub cil_mag: f64,
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub struct CoaParams {
    //TODO: Expand upon LOS system.
    /// Factor controlling to what extent line-of-sight blockage
    /// should be penalized. See SI for further information.
    pub los_penalty: f64,
    /// The distance at which COA signal reaches half-maximum value.
    pub halfmax_dist: f64,
    /// Magnitude of COA that acts on Rac1 activation rates.
    pub vertex_mag: f64,
    //TODO: look up exactly what is being done for this (see where
    // parameter is being generated for hint).
    /// Factor controlling the shape of the exponential modelling
    /// COA interaction (a function shaping parameter). It determines
    /// the distance at which two points would sense COA at half-max
    /// magnitude.
    pub distrib_exp: f64,
    /// If two vertices are within the square root of this distance , then COA cannot occur between
    /// them.
    pub too_close_dist_sq: f64,
}

#[derive(Clone, Copy, Deserialize, Serialize, PartialEq, Debug)]
pub struct ChemAttrParams {
    /// Location of the chemoattractant center.
    pub center: V2d,
    /// Magnitude of chemoattractant a cell would sense if it were
    /// right on top of the chemoattractant source.
    pub center_mag: f64,
    /// Assuming shallow chemoattractant gradient, which can be
    /// modelled using a linear function with slope `slope`.
    pub slope: f64,
}

#[derive(Clone, Deserialize, Serialize, PartialEq, Debug)]
pub struct BdryParams {
    /// Shape of the boundary.
    pub shape: Vec<V2d>,
    /// Bounding box of the boundary.
    pub bbox: BBox,
    /// Should boundary bounding box be checked to see if cell is
    /// within the boundary?
    pub skip_bb_check: bool,
    /// Magnitude of CIL-type interaction.
    pub mag: f64,
}

#[derive(
    Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct InteractionParams {
    pub phys_contact: PhysicalContactParams,
    pub coa: Option<CoaParams>,
    pub chem_attr: Option<ChemAttrParams>,
    pub bdry: Option<BdryParams>,
}

#[derive(
    Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct WorldParameters {
    /// Viscosity value used to calculate change in position of a
    /// vertex due to calculated forces on it.
    pub vertex_eta: f64,
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
#[derive(Clone, Copy, Modify)]
pub struct RawParameters {
    /// Cell diameter.
    pub cell_diam: Length,
    /// Fraction of max force achieved at `rgtp_act_at_max_f`.
    pub halfmax_rgtp_max_f_frac: f64,
    /// Stiffness of the membrane-cortex complex.
    pub stiffness_cortex: Stress,
    /// Typical lamellipod height: typical height of lamellipod (on the order of 100 nm).
    pub lm_h: Length,
    /// Halfmax Rho GTPase activity.
    pub halfmax_rgtp_frac: f64,
    /// Lamellipod stall stress: how much stress can lamellipod exert at most.
    pub lm_ss: Stress,
    /// Friction force opposing RhoA pulling.
    pub rho_friction: f64,
    /// Stiffness of cytoplasm.
    pub stiffness_cyto: Force,
    /// Diffusion rate of Rho GTPase on membrane.
    pub diffusion_rgtp: Diffusion,
    /// Initial distribution of Rac1.
    pub init_rac: RgtpDistribution,
    /// Initial distribution of RhoA.
    pub init_rho: RgtpDistribution,
    /// Baseline Rac1 activation rate.
    pub kgtp_rac: Tinv,
    /// Rac1 auto-activation rate.
    pub kgtp_rac_auto: Tinv,
    /// Baseline Rac1 inactivation rate.
    pub kdgtp_rac: Tinv,
    /// RhoA mediated inhibition of Rac1.
    pub kdgtp_rho_on_rac: Tinv,
    /// Strain at which Rac1 tension-mediated inhibition is half-strength.
    pub halfmax_tension_inhib: f64,
    /// Maximum tension-mediated Rac1 inhibition as a multiple of baseline Rac1 inactivation rate.
    pub tension_inhib: f64,
    /// Rate at which inactive membrane bound Rho GTPase dissociates from the
    /// membrane.
    pub k_mem_off: Tinv,
    /// Rate at which cytosolic Rho GTPase associates with the membrane.
    pub k_mem_on: Tinv,
    /// Baseline RhoA activation rate.
    pub kgtp_rho: Tinv,
    /// RhoA auto-activation rate.
    pub kgtp_auto_rho: Tinv,
    /// Baseline RhoA inactivation rate.
    pub kdgtp_rho: Tinv,
    /// Rac1 mediated inhibition of RhoA.
    pub kdgtp_rac_on_rho: Tinv,
    /// Enable randomization of bursts in Rac1 activity?
    pub randomization: bool,
    /// Average period between randomization events.
    pub rand_avg_t: Time,
    /// Standard deviation of period between randomization events.
    pub rand_std_t: Time,
    /// Magnitude of randomly applied factor affecting Rac1 activation rate: how big a burst?
    pub rand_mag: f64,
    /// Fraction of vertices to be selected for increased Rac1 activation due to random events.
    pub rand_vs: f64,
}

#[derive(
    Copy, Clone, Deserialize, Serialize, Default, Debug, PartialEq,
)]
pub struct Parameters {
    /// Resting cell radius.
    pub cell_r: f64,
    /// Resting edge length.
    pub rest_edge_len: f64,
    /// Resting area.
    pub rest_area: f64,
    /// Stiffness of edge.
    pub stiffness_edge: f64,
    /// Rac1 mediated protrusive force constant.
    pub const_protrusive: f64,
    /// RhoA mediated protrusive force constant.
    pub const_retractive: f64,
    /// Stiffness of cytoplasm.
    pub stiffness_cyto: f64,
    /// Rate of Rho GTPase GDI unbinding and subsequent membrane attachment.
    pub k_mem_on_vertex: f64,
    /// Rate of Rho GTPase membrane disassociation.
    pub k_mem_off: f64,
    /// Diffusion rate of Rho GTPase on membrane.
    pub diffusion_rgtp: f64,
    /// Initial distribution of Rac1.
    pub init_rac: RgtpDistribution,
    /// Initial distribution of RhoA.
    pub init_rho: RgtpDistribution,
    /// Halfmax Rho GTPase activity per vertex.
    pub halfmax_vertex_rgtp: f64,
    /// Halfmax Rho GTPase activity per vertex as concentration.
    pub halfmax_vertex_rgtp_conc: f64,
    /// Baseline Rac1 activation rate.
    pub kgtp_rac: f64,
    /// Rac1 auto-activation rate as a multiple of baseline Rac1 activation rate.
    pub kgtp_rac_auto: f64,
    /// Baseline Rac1 inactivation rate.
    pub kdgtp_rac: f64,
    /// RhoA mediated inhibition of Rac1 as a multiple of baseline Rac1 inactivation rate.
    pub kdgtp_rho_on_rac: f64,
    /// Strain at which Rac1 tension-mediated inhibition is half-strength.
    pub halfmax_tension_inhib: f64,
    /// Tension-mediated Rac1 inhibition as a multiple of baseline Rac1 inactivation rate.
    pub tension_inhib: f64,
    /// Baseline RhoA activation rate.
    pub kgtp_rho: f64,
    /// RhoA auto-activation rate as a multiple of baseline RhoA activation rate.
    pub kgtp_rho_auto: f64,
    /// Baseline RhoA inactivation rate.
    pub kdgtp_rho: f64,
    /// Rac1 mediated inhibition of RhoA as a multiple of baseline RhoA inactivation rate.
    pub kdgtp_rac_on_rho: f64,
    /// Enable randomization of bursts in Rac1 activity?
    pub randomization: bool,
    /// Average time between random events, in timepoints.
    pub rand_avg_t: f64,
    /// Standard deviation of time between random events, in timepoints.
    pub rand_std_t: f64,
    /// Magnitude of factor randomly applied to Rac1 activation rate.
    pub rand_mag: f64,
    /// Number of vertices to be selected for random Rac1 activity boost.
    pub num_rand_vs: u32,
}

impl RawParameters {
    pub fn refine(&self, bq: &CharQuantities) -> Parameters {
        let cell_r = self.cell_diam.scale(0.5);
        let rel = self.cell_diam.scale((PI / (NVERTS as f64)).sin());
        let ra = Length(1.0)
            .pow(2.0)
            .scale(calc_init_cell_area(cell_r.number()));
        let const_protrusive =
            (self.lm_h.g() * self.lm_ss.g() * rel.g())
                .scale(self.halfmax_rgtp_max_f_frac);
        let const_retractive =
            const_protrusive.scale(self.rho_friction);
        let halfmax_vertex_rgtp =
            self.halfmax_rgtp_frac / NVERTS as f64;
        let halfmax_vertex_rgtp_conc =
            rel.pow(-1.0).scale(halfmax_vertex_rgtp);
        let stiffness_edge = self.stiffness_cortex.g() * bq.l3d.g();
        let stiffness_cyto =
            self.stiffness_cyto.g().scale(1.0 / NVERTS as f64);

        Parameters {
            cell_r: bq.normalize(&cell_r),
            rest_edge_len: bq.normalize(&rel),
            rest_area: bq.normalize(&ra),
            stiffness_edge: bq.normalize(&stiffness_edge),
            const_protrusive: bq.normalize(&const_protrusive),
            const_retractive: bq.normalize(&const_retractive),
            stiffness_cyto: bq.normalize(&stiffness_cyto),
            k_mem_on_vertex: bq.normalize(&self.k_mem_on)
                / NVERTS as f64,
            k_mem_off: bq.normalize(&self.k_mem_off),
            diffusion_rgtp: bq.normalize(&self.diffusion_rgtp),
            init_rac: self.init_rac,
            init_rho: self.init_rho,
            halfmax_vertex_rgtp,
            halfmax_vertex_rgtp_conc: bq
                .normalize(&halfmax_vertex_rgtp_conc),
            kgtp_rac: bq.normalize(&self.kgtp_rac),
            kgtp_rac_auto: bq.normalize(&self.kgtp_rac_auto),
            kdgtp_rac: bq.normalize(&self.kdgtp_rac),
            kdgtp_rho_on_rac: bq.normalize(&self.kdgtp_rho_on_rac),
            halfmax_tension_inhib: self.halfmax_tension_inhib,
            tension_inhib: self.tension_inhib,
            kgtp_rho: bq.normalize(&self.kgtp_rho),
            kgtp_rho_auto: bq.normalize(&self.kgtp_auto_rho),
            kdgtp_rho: bq.normalize(&self.kdgtp_rho),
            kdgtp_rac_on_rho: bq.normalize(&self.kdgtp_rac_on_rho),
            randomization: self.randomization,
            rand_avg_t: bq.normalize(&self.rand_avg_t).ceil(),
            rand_std_t: bq.normalize(&self.rand_std_t).ceil(),
            rand_mag: self.rand_mag,
            num_rand_vs: (self.rand_vs * NVERTS as f64) as u32,
        }
    }
}

/// Calculate the area of an "ideal" initial cell of radius R, if it has
/// `NVERTS`  vertices.
pub fn calc_init_cell_area(r: f64) -> f64 {
    let poly_coords = (0..NVERTS)
        .map(|vix| {
            let theta = (vix as f64) / (NVERTS as f64) * 2.0 * PI;
            V2d {
                x: r * theta.cos(),
                y: r * theta.sin(),
            }
        })
        .collect::<Vec<V2d>>();
    calc_poly_area(&poly_coords)
}
