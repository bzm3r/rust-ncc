use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::experiments::{
    gen_default_char_quants, gen_default_phys_contact_dist,
    gen_default_raw_params, gen_default_vertex_viscosity, CellGroup,
    Experiment, GroupBBox,
};
use crate::math::v2d::V2D;
use crate::parameters::quantity::{Length, Quantity};
use crate::parameters::{
    CharQuantities, RawCloseBounds, RawInteractionParams,
    RawParameters, RawPhysicalContactParams, RawWorldParameters,
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
        gen_raw_params(rng, false, group0_marked, group1_marked);
    let raw_params1 =
        gen_raw_params(rng, false, group1_marked, group0_marked);
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
        (Length(0.0), raw_params1.cell_diam.mul_number(1.0));
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
    char_quants: &CharQuantities,
) -> RawWorldParameters {
    // Some(RawCoaParams {
    //     los_penalty: 2.0,
    //     range: Length(100.0).micro(),
    //     mag: 100.0,
    // })
    let one_at = gen_default_phys_contact_dist();
    RawWorldParameters {
        vertex_eta: gen_default_vertex_viscosity(char_quants),
        interactions: RawInteractionParams {
            coa: None,
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: RawCloseBounds::new(
                    one_at.mul_number(2.0),
                    one_at,
                ),
                adh_mag: None,
                cal_mag: None,
                cil_mag: 60.0,
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
    let char_quants = gen_default_char_quants();
    let world_parameters =
        raw_world_parameters(&char_quants).refine(&char_quants);
    let cell_groups = cell_groups(&mut rng, &char_quants);
    Experiment {
        file_name: "cil_test".to_string(),
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
            frac: 0.0,
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
