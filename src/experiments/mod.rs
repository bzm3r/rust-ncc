// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod pairs;
pub mod single;

use crate::cell::chemistry::{DistributionScheme, DistributionType, RgtpDistribution};
use crate::math::v2d::V2d;
use crate::parameters::quantity::{Force, Length, Quantity, Stress, Time, Tinv, Viscosity};
use crate::parameters::{CharQuantities, Parameters, RawParameters, WorldParameters};
use crate::NVERTS;

/// Specifies initial placement of the group.
pub struct GroupLayout {
    /// Width of group in terms of number of cells.
    pub width: u32,
    /// Height of group in terms of number if cells.
    pub height: u32,
    /// Bottom left of the group in micrometers.
    pub bottom_left: V2d,
}

/// Information required for a cell group to be created.
pub struct CellGroup {
    /// The number of cells in the group.
    pub num_cells: u32,
    /// Initial layout of the cell group.
    pub layout: GroupLayout,
    /// Parameters shared by all cells in this group.
    pub parameters: Parameters,
}

/// Information required to create an experiment.
pub struct Experiment {
    pub title: String,
    /// Characteristic quantities.
    pub char_quants: CharQuantities,
    pub world_parameters: WorldParameters,
    /// List of cell groups involved in this experiment.
    pub cell_groups: Vec<CellGroup>,
}

/// Generate default characteristic quantities.
/// If you need to use different characteristic quantities
/// than the ones provided in the defaults, implement a
/// function that generates `CharQuantities` within your
/// particular experiment file. You can use this function
/// as a template.
///
/// Refer to SI of first two papers for justification of
/// the values used.
// TODO: Document all the justifications here, rather than having to refer
// to the SI.
fn gen_default_char_quants() -> CharQuantities {
    // Stress on lamellipod is on order of 1kPa, height of lamellipod on order of 100 nm, length of edge on order of 10 um
    let f = (Stress(1.0).kilo().g() * Length(100.0).nano().g() * Length(10.0).micro().g())
        .to_force()
        .unwrap();

    CharQuantities {
        eta: Viscosity(0.1),
        l: Length(1.0).micro(),
        t: Time(2.0),
        f,
        l3d: Length(10e-6),
        k_mem_off: Tinv(0.15),
        k_mem_on: Tinv(0.02),
        kgtp: Tinv(1e-4),
        kdgtp: Tinv(1e-4),
        frac_rgtp: 0.1,
    }
}

/// Generate default raw parameters for cells.
/// If you need to use different raw parameters
/// than the ones provided in the defaults, implement a
/// function that generates `RawParameters` within your
/// particular experiment file. You can use this function
/// as a template.
fn gen_default_raw_params() -> RawParameters {
    let rgtp_d = (Length(0.1_f32.sqrt()).micro().pow(2.0).g() / Time(1.0).g())
        .to_diffusion()
        .unwrap();
    let init_rac = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
    )
    .unwrap();
    RawParameters {
        cell_diam: Length(40.0).micro(),
        stiffness_cortex: Stress(8.0).kilo(),
        lm_h: Length(200.0).nano(),
        halfmax_rgtp_max_f_frac: 0.3,
        halfmax_rgtp_frac: 0.4,
        lm_ss: Stress(10.0).kilo(),
        rho_friction: 0.2,
        stiffness_ctyo: Force(1e-7),
        diffusion_rgtp: rgtp_d,
        tot_rac: 2.5e6,
        tot_rho: 1e6,
        kgtp_rac: 24.0,
        kgtp_rac_auto: 500.0,
        kdgtp_rac: 8.0,
        kdgtp_rho_on_rac: 4000.0,
        halfmax_tension_inhib: 0.1,
        tension_inhib: 40.0,
        kgtp_rho: 28.0,
        kgtp_auto_rho: 390.0,
        kdgtp_rho: 60.0,
        kdgtp_rac_on_rho: 400.0,
        randomization: true,
        rand_avg_t: Time(40.0 * 60.0),
        rand_std_t: Time(0.2 * 40.0 * 60.0),
        rand_mag: 10.0,
        rand_vs: 0.25,
        init_rac,
        init_rho: init_rac,
    }
}

/// We take the viscosity of the world to be 0.29 N m^-2. We
/// divide viscosity by the number of vertices, on a cell in
/// order to scale it properly.
///
/// See SI for justification.
// TODO: put justification here.
fn gen_default_viscosity() -> Viscosity {
    Viscosity(0.29).mul_number(1.0 / (NVERTS as f32))
}
