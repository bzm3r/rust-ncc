#![allow(unused)]
use crate::cell::chemistry::{DistributionScheme, DistributionType, RgtpDistribution};
use crate::experiments::{CellGroup, Experiment, GroupLayout};
use crate::interactions::CilMat;
use crate::math::p2d::P2D;
use crate::parameters::quantity::{Force, Length, Quantity, Stress, Time, Tinv, Viscosity};
use crate::parameters::{BasicQuants, RawParameters, RawWorldParameters};
use crate::NVERTS;

fn group_layout(bq: &BasicQuants) -> GroupLayout {
    let raw_centroid = [Length(0.0), Length(0.0)];
    let centroid = P2D {
        x: bq.normalize(&raw_centroid[0]),
        y: bq.normalize(&raw_centroid[1]),
    };
    GroupLayout {
        width: 2,
        height: 1,
        bottom_left: centroid,
    }
}

fn cell_groups(bq: &BasicQuants) -> Vec<CellGroup> {
    vec![CellGroup {
        num_cells: 1,
        layout: group_layout(bq),
        parameters: raw_parameters().gen_parameters(bq),
    }]
}

fn basic_quants() -> BasicQuants {
    // Stress on lamellipod is on order of 1kPa, height of lamellipod on order of 100 nm, length of edge on order of 10 um
    let f = (Stress(1.0).kilo().g() * Length(100.0).nano().g() * Length(10.0).micro().g())
        .to_force()
        .unwrap();
    let eta = Viscosity(0.1);

    BasicQuants {
        eta,
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

fn gen_cil_mat() -> CilMat {
    CilMat::new(2, 60.0)
}

fn raw_world_parameters() -> RawWorldParameters {
    RawWorldParameters {
        vertex_eta: Viscosity(0.29).mulf(1.0 / (NVERTS as f32)),
        close_criterion: Length(0.5).micro(),
        cil: gen_cil_mat(),
    }
}

fn raw_parameters() -> RawParameters {
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
        randomization: true,
        rand_avg_t: Time(40.0 * 60.0),
        rand_std_t: Time(0.2 * 40.0 * 60.0),
        rand_mag: 10.0,
        rand_vs: 0.25,
        init_rac,
        init_rho: init_rac,
    }
}

pub fn generate() -> Experiment {
    let basic_quants = basic_quants();
    let world_parameters = raw_world_parameters().normalize(&basic_quants);
    let cell_groups = cell_groups(&basic_quants);
    Experiment {
        title: "CIL test".to_string(),
        basic_quants,
        world_parameters,
        cell_groups,
    }
}
