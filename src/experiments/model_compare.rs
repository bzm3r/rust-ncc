use crate::experiments::{
    gen_default_char_quants, gen_default_phys_contact_dist,
    gen_default_raw_params, gen_default_vertex_viscosity, CellGroup,
    Experiment, GroupBBox,
};
use crate::math::v2d::V2D;
use crate::parameters::quantity::{Length, Quantity};
use crate::parameters::{
    CharQuantities, RawCloseBounds, RawCoaParams,
    RawInteractionParams, RawParameters, RawPhysicalContactParams,
    RawWorldParameters,
};
use crate::utils::pcg32::Pcg32;

use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::NVERTS;
use rand::SeedableRng;

/// Generate the group layout to use for this experiment.
fn group_layout(
    group_ix: usize,
    char_quants: &CharQuantities,
) -> Result<GroupBBox, String> {
    // specify initial location of group centroid
    let centroid = V2D {
        x: char_quants
            .normalize(&Length((group_ix as f64) * 40.0).micro().g()),
        y: char_quants.normalize(&Length(0.0)),
    };
    let r = GroupBBox {
        width: 1,
        height: 1,
        bottom_left: centroid,
    };
    if r.width * r.height < 1 {
        Err(String::from(
            "Group layout area is too small to contain required number of cells.",
        ))
    } else {
        Ok(r)
    }
}

fn gen_raw_params(
    group_ix: usize,
    rng: &mut Pcg32,
    randomization: bool,
) -> RawParameters {
    let mut right = [false; NVERTS];
    right.iter_mut().enumerate().for_each(|(i, x)| match i {
        0 | 1 | 2 | 3 => {
            *x = true;
        }
        _ => {}
    });
    let mut left = [false; NVERTS];
    left.iter_mut().enumerate().for_each(|(i, x)| match i {
        8 | 9 | 10 | 11 => {
            *x = true;
        }
        _ => {}
    });
    let (specific_rac, specific_rho) = match group_ix {
        0 => (right, left),
        1 => (left, right),
        _ => panic!("received group ix > 1"),
    };
    let init_rac = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.3,
            ty: DistributionType::SpecificUniform(specific_rac),
        },
        DistributionScheme {
            frac: 0.3,
            ty: DistributionType::SpecificUniform(specific_rac),
        },
        rng,
    )
    .unwrap();
    let init_rho = RgtpDistribution::generate(
        DistributionScheme {
            frac: 0.3,
            ty: DistributionType::SpecificUniform(specific_rho),
        },
        DistributionScheme {
            frac: 0.3,
            ty: DistributionType::SpecificUniform(specific_rho),
        },
        rng,
    )
    .unwrap();

    let mut raw_params = gen_default_raw_params(rng, randomization);
    raw_params.modify_init_rac(init_rac);
    raw_params.modify_init_rho(init_rho);

    raw_params
}

/// Define the cell groups that will exist in this experiment.
fn cell_groups(
    rng: &mut Pcg32,
    randomization: bool,
    cq: &CharQuantities,
    num_cells_per_group: Vec<usize>,
) -> Vec<CellGroup> {
    num_cells_per_group
        .iter()
        .enumerate()
        .filter_map(|(group_ix, &num_cells)| {
            if num_cells > 0 {
                Some(CellGroup {
                    num_cells,
                    layout: group_layout(group_ix, cq).unwrap(),
                    parameters: gen_raw_params(
                        group_ix,
                        rng,
                        randomization,
                    )
                    .gen_parameters(cq),
                })
            } else {
                None
            }
        })
        .collect()
}

/// Generate raw world parameters, in particular, how
/// cells interact with each other, and any boundaries.
fn raw_world_parameters(
    char_quants: &CharQuantities,
    coa_mag: Option<f64>,
    cil_mag: f64,
) -> RawWorldParameters {
    // Some(RawCoaParams {
    //     los_penalty: 2.0,
    //     range: Length(100.0).micro(),
    //     mag: 100.0,
    // })
    let one_at = gen_default_phys_contact_dist();
    // Some(RawCoaParams {
    //                 los_penalty: 2.0,
    //                 range: Length(220.0).micro(),
    //                 mag: 24.0,
    //             })
    let _vertex_eta = gen_default_vertex_viscosity(char_quants);
    RawWorldParameters {
        vertex_eta: gen_default_vertex_viscosity(char_quants),
        interactions: RawInteractionParams {
            coa: coa_mag.map(|x| RawCoaParams {
                los_penalty: 2.0,
                halfmax_dist: Length(110.0).micro(),
                mag: x,
            }),
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: RawCloseBounds::new(
                    one_at.mul_number(3.0),
                    one_at,
                ),
                adh_mag: None,
                cal_mag: None,
                cil_mag,
            },
        },
    }
}

/// Generate the experiment, so that it can be run.
pub fn generate(
    seed: Option<u64>,
    num_cells_per_group: Vec<usize>,
    randomization: bool,
    coa_mag: Option<f64>,
    cil_mag: f64,
) -> Experiment {
    let mut rng = match seed {
        Some(s) => Pcg32::seed_from_u64(s),
        None => Pcg32::from_entropy(),
    };
    let char_quants = gen_default_char_quants();
    let world_parameters =
        raw_world_parameters(&char_quants, coa_mag, cil_mag)
            .refine(&char_quants);
    let cell_groups = cell_groups(
        &mut rng,
        randomization,
        &char_quants,
        num_cells_per_group,
    );
    Experiment {
        file_name: "n_cells".to_string(),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}
