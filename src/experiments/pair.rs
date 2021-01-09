#![allow(unused)]
use rand::SeedableRng;
use rust_ncc::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use rust_ncc::experiments::{
    gen_default_adhesion_mag, gen_default_char_quants,
    gen_default_phys_contact_dist, gen_default_raw_params,
    gen_default_viscosity, CellGroup, Experiment, GroupBBox,
};
use rust_ncc::interactions::dat_sym2d::SymCcDat;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::quantity::{Force, Length, Quantity};
use rust_ncc::parameters::{
    CharQuantities, CoaParams, PhysicalContactParams, RawCloseBounds,
    RawCoaParams, RawInteractionParams, RawParameters,
    RawPhysicalContactParams, RawWorldParameters,
};
use rust_ncc::utils::pcg32::Pcg32;
use rust_ncc::NVERTS;

/// Generate the group layout to use for this experiment.
fn group_layout(
    num_cells: usize,
    char_quants: &CharQuantities,
) -> Result<GroupBBox, String> {
    // specify initial location of group centroid
    let centroid = V2D {
        x: char_quants.normalize(&Length(0.0)),
        y: char_quants.normalize(&Length(0.0)),
    };
    let r = GroupBBox {
        width: 1,
        height: 2,
        bottom_left: centroid,
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
    let num_cells = 2;
    vec![CellGroup {
        num_cells,
        layout: group_layout(num_cells, cq).unwrap(),
        parameters: gen_default_raw_params(rng, true)
            .gen_parameters(cq),
    }]
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
    char_quants: &CharQuantities,
) -> RawWorldParameters {
    // Some(RawCoaParams {
    //     los_penalty: 2.0,
    //     range: Length(100.0).micro(),
    //     mag: 100.0,
    // })
    let one_at = gen_default_phys_contact_dist();
    RawWorldParameters {
        vertex_eta: gen_default_viscosity(),
        interactions: RawInteractionParams {
            coa: None,
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: RawCloseBounds::new(
                    one_at.mul_number(2.0),
                    one_at,
                ),
                adh_mag: Some(gen_default_adhesion_mag(
                    char_quants,
                    0.0,
                )),
                cal_mag: Some(60.0),
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
        file_name: "pair".to_string(),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}
