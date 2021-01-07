#![allow(unused)]
use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::experiments::{
    gen_default_char_quants, gen_default_phys_contact_dist,
    gen_default_raw_params, gen_default_viscosity, CellGroup,
    Experiment, GroupBBox,
};
use crate::interactions::dat_sym2d::SymCcDat;
use crate::math::v2d::V2D;
use crate::parameters::quantity::Length;
use crate::parameters::{
    CharQuantities, RawInteractionParams, RawParameters,
    RawPhysicalContactParams, RawWorldParameters,
};
use crate::utils::pcg32::Pcg32;
use crate::NVERTS;
use rand::SeedableRng;

/// Generate the group layout to use for this experiment.
fn group_layout(
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
    let num_cells = 1;
    vec![CellGroup {
        num_cells,
        layout: group_layout(num_cells, cq).unwrap(),
        parameters: gen_default_raw_params(rng, true)
            .gen_parameters(cq),
    }]
}

/// Generate CAL values between different cells.
fn gen_cal_mat() -> SymCcDat<f64> {
    SymCcDat::<f64>::new(2, 60.0)
}

/// Generate CIL values between different cells.
fn gen_cil_mat() -> SymCcDat<f64> {
    SymCcDat::<f64>::new(2, 60.0)
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
