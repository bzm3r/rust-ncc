#![allow(unused)]
use rand::SeedableRng;
use rust_ncc::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use rust_ncc::experiments::{
    gen_default_char_quants, gen_default_phys_contact_dist,
    gen_default_raw_params, gen_default_viscosity, CellGroup,
    Experiment, GroupBBox,
};
use rust_ncc::interactions::dat_sym2d::SymCcDat;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::quantity::Length;
use rust_ncc::parameters::{
    CharQuantities, RawInteractionParams, RawParameters,
    RawPhysicalContactParams, RawWorldParameters,
};
use rust_ncc::utils::pcg32::Pcg32;
use rust_ncc::NVERTS;

/// Generate the group layout to use for this experiment.
fn group_layout(
    num_cells: u32,
    char_quants: &CharQuantities,
    bottom_left: (Length, Length),
    width: u32,
    height: u32,
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
    let num_cells = 1;
    vec![CellGroup {
        num_cells,
        layout: group_layout(num_cells, cq).unwrap(),
        parameters: gen_default_raw_params(rng, true)
            .gen_parameters(cq),
    }]
}

/// Generate CAL values between different cells.
fn gen_cal_mat() -> SymCcDat<f32> {
    SymCcDat::<f32>::new(2, 60.0)
}

/// Generate CIL values between different cells.
fn gen_cil_mat() -> SymCcDat<f32> {
    SymCcDat::<f32>::new(2, 60.0)
}

/// Generate raw world parameters, in particular, how
/// cells interact with each other, and any boundaries.
fn raw_world_parameters() -> RawWorldParameters {
    RawWorldParameters {
        vertex_eta: gen_default_viscosity(),
        interactions: RawInteractionParams {
            coa: None,
            chem_attr: None,
            bdry: None,
            phys_contact: RawPhysicalContactParams {
                range: gen_default_phys_contact_dist(),
                adh_mag: None,
                cal_mag: None,
                cil_mag: 0.0,
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
        raw_world_parameters().refine(&char_quants);
    let cell_groups = cell_groups(&mut rng, &char_quants);
    Experiment {
        file_name: "single".to_string(),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}
