#![allow(unused)]
use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::experiments::{
    gen_default_adhesion_mag, gen_default_char_quants,
    gen_default_phys_contact_dist, gen_default_viscosity, CellGroup,
    Experiment, GroupBBox,
};
use crate::interactions::dat_sym2d::SymCcDat;
use crate::math::v2d::V2D;
use crate::parameters::quantity::{
    Force, Length, Quantity, Stress, Time,
};
use crate::parameters::{
    CharQuantities, CoaParams, PhysicalContactParams, RawCloseBounds,
    RawCoaParams, RawInteractionParams, RawParameters,
    RawPhysicalContactParams, RawWorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::NVERTS;
use rand::SeedableRng;

/// Generate the group layout to use for this experiment.
fn group_bbox(
    num_cells: usize,
    char_quants: &CharQuantities,
    bottom_left: (Length, Length),
    width: usize,
    height: usize,
) -> Result<GroupBBox, String> {
    // specify initial location of group bottom left
    let bottom_left = V2D {
        x: char_quants.normalize(&bottom_left.0),
        y: char_quants.normalize(&bottom_left.1),
    };
    let r = GroupBBox {
        width,
        height,
        bottom_left,
    };

    if r.width * r.height > num_cells {
        Err(String::from(
            "Group layout area is too small to contain required number of cells.",
        ))
    } else {
        Ok(r)
    }
}

/// Define the cell groups that will exist in this experiment.
fn cell_groups(
    rng: &mut Pcg32,
    cq: &CharQuantities,
) -> Vec<CellGroup> {
    let group0_marked_rho = [
        false, true, true, true, true, true, true, true, false,
        false, false, false, false, false, false, false,
    ];
    let group0_marked = [
        false, false, false, false, false, false, false, false,
        false, true, true, true, true, true, true, true,
    ];
    let group1_marked = [
        false, false, false, false, false, false, false, false,
        false, true, true, true, true, true, true, true,
    ];
    let raw_params0 = gen_default_raw_params_0(
        rng,
        false,
        group0_marked,
        group0_marked_rho,
    );
    let raw_params1 = gen_default_raw_params_1(
        rng,
        false,
        group1_marked,
        group0_marked,
    );
    let params0 = raw_params0.gen_parameters(cq);
    let params1 = raw_params1.gen_parameters(cq);
    let bottom_left0 = (Length(0.0), Length(0.0));
    let num_cells0 = 1;
    let group0_layout = CellGroup {
        num_cells: num_cells0,
        layout: group_bbox(num_cells0, cq, bottom_left0, 1, 1)
            .unwrap(),
        parameters: params0,
    };
    let bottom_left1 =
        (Length(0.0), raw_params1.cell_diam.mul_number(1.1));
    let num_cells1 = 1;
    let group1_layout = CellGroup {
        num_cells: num_cells1,
        layout: group_bbox(num_cells1, cq, bottom_left1, 1, 1)
            .unwrap(),
        parameters: params1,
    };
    vec![group0_layout, group1_layout]
}

/// Generate CAL values between different cells.
fn gen_cal_mat() -> SymCcDat<f32> {
    SymCcDat::<f32>::new(2, 0.0)
}

/// Generate CIL values between different cells (see SI for
/// justification).
fn gen_cil_mat() -> SymCcDat<f32> {
    SymCcDat::<f32>::new(2, 60.0)
}

/// Generate raw world parameters, in particular, how
/// cells interact with each other, and any boundaries.
fn raw_world_parameters(
    coa_mag: Option<f32>,
    char_quants: &CharQuantities,
) -> RawWorldParameters {
    // Some(RawCoaParams {
    //     los_penalty: 2.0,
    //     range: Length(100.0).micro(),
    //     mag: 100.0,
    // })
    let one_at = gen_default_phys_contact_dist();
    let coa = RawCoaParams::default_with_mag(coa_mag);
    RawWorldParameters {
        vertex_eta: gen_default_viscosity(),
        interactions: RawInteractionParams {
            coa,
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: RawCloseBounds::new(
                    one_at.mul_number(2.0),
                    one_at,
                ),
                // adh_mag: None,
                adh_mag: Some(gen_default_adhesion_mag(
                    char_quants,
                    10.0,
                )),
                cal_mag: Some(10.0),
                cil_mag: 10.0,
            },
        },
    }
}

/// Generate the experiment, so that it can be run.
pub fn generate(seed: Option<u64>) -> Experiment {
    let mut rng = match seed {
        Some(s) => Pcg32::seed_from_u64(s),
        None => Pcg32::from_entropy(),
    };
    let cil = 10.0;
    let cal: Option<f32> = Some(10.0);
    let adh = Some(10.0);
    let coa: Option<f32> = None;

    let char_quants = gen_default_char_quants();
    let world_parameters =
        raw_world_parameters(None, &char_quants).refine(&char_quants);
    let cell_groups = cell_groups(&mut rng, &char_quants);

    //convert the option into string
    let cal = if let Some(i) = cal {
        i.to_string()
    } else {
        "None".to_string()
    };

    let adh = if let Some(i) = adh {
        i.to_string()
    } else {
        "None".to_string()
    };

    let coa = if let Some(i) = coa {
        i.to_string()
    } else {
        "None".to_string()
    };

    Experiment {
        // file_name: format!(
        //     "cal_test_cil={}_cal={}_adh={}",
        //     cil, cal, adh
        // ),
        file_name: format!(
            "cal_test_cil={}_cal={}_adh={}_coa={}",
            cil, cal, adh, coa
        ),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}

fn gen_default_raw_params_0(
    rng: &mut Pcg32,
    randomization: bool,
    marked_rac: [bool; NVERTS],
    marked_rho: [bool; NVERTS],
) -> RawParameters {
    //println!("marking: {:?}", &marked_rac);

    let rgtp_d = (Length(0.1_f32.sqrt()).micro().pow(2.0).g()
        / Time(1.0).g())
    .to_diffusion()
    .unwrap();

    let init_rac = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.03,
            ty: DistributionType::SpecificUniform(marked_rac),
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        rng,
    )
    .unwrap();

    let init_rho = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::SpecificUniform(marked_rho),
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        rng,
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
        randomization,
        rand_avg_t: Time(40.0 * 60.0),
        rand_std_t: Time(0.2 * 40.0 * 60.0),
        rand_mag: 10.0,
        rand_vs: 0.25,
        init_rac,
        init_rho,
    }
}

fn gen_default_raw_params_1(
    rng: &mut Pcg32,
    randomization: bool,
    marked_rac: [bool; NVERTS],
    marked_rho: [bool; NVERTS],
) -> RawParameters {
    //println!("marking: {:?}", &marked_rac);

    let rgtp_d = (Length(0.1_f32.sqrt()).micro().pow(2.0).g()
        / Time(1.0).g())
    .to_diffusion()
    .unwrap();

    let init_rac = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::SpecificUniform(marked_rac),
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        rng,
    )
    .unwrap();

    let init_rho = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::SpecificUniform(marked_rho),
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        rng,
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
        randomization,
        rand_avg_t: Time(40.0 * 60.0),
        rand_std_t: Time(0.2 * 40.0 * 60.0),
        rand_mag: 10.0,
        rand_vs: 0.25,
        init_rac,
        init_rho,
    }
}
