#![allow(unused)]
use crate::cell::chemistry::{
    DistributionScheme, DistributionType, RgtpDistribution,
};
use crate::experiments::{
    gen_default_adhesion_mag, gen_default_char_quants,
    gen_default_phys_contact_dist, gen_default_raw_params,
    gen_default_viscosity, CellGroup, Experiment, GroupLayout,
};
use crate::interactions::symdat2d::SymCcDat;
use crate::math::v2d::V2d;
use crate::parameters::quantity::{Force, Length};
use crate::parameters::{
    CharQuantities, PhysicalContactParams, RawInteractionParams,
    RawParameters, RawPhysicalContactParams, RawWorldParameters,
};
use crate::NVERTS;
use rand::SeedableRng;
use rand_pcg::Pcg64;

/// Generate the group layout to use for this experiment.
fn group_layout(
    num_cells: u32,
    char_quants: &CharQuantities,
) -> Result<GroupLayout, String> {
    // specify initial location of group centroid
    let centroid = V2d {
        x: char_quants.normalize(&Length(0.0)),
        y: char_quants.normalize(&Length(0.0)),
    };
    let r = GroupLayout {
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
    rng: &mut Pcg64,
    cq: &CharQuantities,
) -> Vec<CellGroup> {
    let num_cells = 2;
    vec![CellGroup {
        num_cells,
        layout: group_layout(num_cells, cq).unwrap(),
        parameters: gen_default_raw_params(rng, false)
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
fn raw_world_parameters() -> RawWorldParameters {
    RawWorldParameters {
        vertex_eta: gen_default_viscosity(),
        interactions: RawInteractionParams {
            coa: None,
            chem_attr: None,
            bdry: None,
            phys_contact: Some(RawPhysicalContactParams {
                range: gen_default_phys_contact_dist(),
                adh_mag: gen_default_adhesion_mag(),
                cal_mag: 0.0,
                cil_mag: 60.0,
            }),
        },
    }
}

/// Generate the experiment, so that it can be run.
pub fn generate(seed: Option<u64>) -> Experiment {
    let mut rng = match seed {
        Some(s) => Pcg64::seed_from_u64(s),
        None => Pcg64::from_entropy(),
    };
    let char_quants = gen_default_char_quants();
    let world_parameters =
        raw_world_parameters().refine(&char_quants);
    let cell_groups = cell_groups(&mut rng, &char_quants);
    Experiment {
        title: "a pair of cells".to_string(),
        char_quants,
        world_parameters,
        cell_groups,
        rng,
        seed,
    }
}
