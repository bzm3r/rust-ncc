use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::experiments::{
    gen_default_adhesion_mag, gen_default_char_quants,
    gen_default_phys_contact_dist, gen_default_raw_params,
    gen_default_vertex_viscosity, CellGroup, Experiment, GroupBBox,
};
use crate::math::v2d::V2D;
use crate::parameters::quantity::{Length, Quantity};
use crate::parameters::{
    CharQuantities, RawCloseBounds, RawCoaParams,
    RawInteractionParams, RawParameters, RawPhysicalContactParams,
    RawWorldParameters,
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
    let group0_marked = [
        false, true, true, true, true, true, true, true, false,
        false, false, false, false, false, false, false,
    ];
    let group1_marked = [
        false, false, false, false, false, false, false, false,
        false, true, true, true, true, true, true, true,
    ];
    let raw_params0 =
        gen_raw_params(rng, true, group0_marked, group1_marked);
    let raw_params1 =
        gen_raw_params(rng, true, group1_marked, group0_marked);
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
        (Length(0.0), raw_params1.cell_diam.mul_number(2.0));
    let num_cells1 = 1;
    let group1_layout = CellGroup {
        num_cells: num_cells1,
        layout: group_bbox(num_cells1, cq, bottom_left1, 1, 1)
            .unwrap(),
        parameters: params1,
    };
    vec![group0_layout, group1_layout]
}

/// Generate raw world parameters, in particular, how
/// cells interact with each other, and any boundaries.
fn raw_world_parameters(
    coa_mag: Option<f64>,
    adh_mag: Option<f64>,
    cal_mag: Option<f64>,
    cil_mag: f64,
    char_quants: &CharQuantities,
) -> RawWorldParameters {
    // Some(RawCoaParams {
    //     los_penalty: 2.0,
    //     range: Length(100.0).micro(),
    //     mag: 100.0,
    // })
    let one_at = gen_default_phys_contact_dist();
    let coa = RawCoaParams::default_with_mag(coa_mag);
    let adh_mag = if let Some(x) = adh_mag {
        Some(gen_default_adhesion_mag(char_quants, x))
    } else {
        None
    };
    RawWorldParameters {
        vertex_eta: gen_default_vertex_viscosity(char_quants),
        interactions: RawInteractionParams {
            coa,
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: RawCloseBounds::new(
                    one_at.mul_number(2.0),
                    one_at,
                ),
                adh_mag,
                cal_mag,
                cil_mag,
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
    let cil = 60.0;
    let cal: Option<f64> = None;
    let adh: Option<f64> = Some(10.0);
    let coa: Option<f64> = Some(24.0);

    let char_quants = gen_default_char_quants();
    let world_parameters =
        raw_world_parameters(coa, adh, cal, cil, &char_quants)
            .refine(&char_quants);
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

    let seed_string = if let Some(i) = seed {
        i.to_string()
    } else {
        "None".to_string()
    };

    Experiment {
        file_name: format!(
            "separated_pair_cil={}_cal={}_adh={}_coa={}_seed={}",
            cil, cal, adh, coa, seed_string
        ),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}

fn gen_raw_params(
    rng: &mut Pcg32,
    randomization: bool,
    marked_rac: [bool; NVERTS],
    marked_rho: [bool; NVERTS],
) -> RawParameters {
    let init_rac = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::SpecificRandom(marked_rac),
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
            ty: DistributionType::SpecificRandom(marked_rho),
        },
        DistributionScheme {
            frac: 0.1,
            ty: DistributionType::Random,
        },
        rng,
    )
    .unwrap();

    let mut raw_params = gen_default_raw_params(rng, randomization);
    raw_params.modify_init_rac(init_rac);
    raw_params.modify_init_rho(init_rho);

    raw_params
}
