use crate::cell::chemistry::{distrib_gens, RgtpDistribution};
use crate::exp_setup::defaults::RAW_COA_PARAMS_WITH_ZERO_MAG;
use crate::exp_setup::exp_parser::ExperimentArgs;
use crate::exp_setup::markers::mark_between_angles;
use crate::exp_setup::{
    defaults, CellGroup, Experiment, ExperimentType, GroupBBox,
};
use crate::math::radians::{Radians, RAD_0, RAD_PI};
use crate::math::v2d::V2d;
use crate::parameters::quantity::{Length, Quantity};
use crate::parameters::{
    CharQuantities, RawCloseBounds, RawInteractionParams,
    RawParameters, RawPhysicalContactParams,
};
use crate::utils::pcg32::Pcg32;
use crate::Directories;
use rand::SeedableRng;

/// Generate the group bounding box to use for this experiment.
fn group_bbox(
    num_cells: usize,
    char_quants: &CharQuantities,
    bottom_left: (Length, Length),
    width: usize,
    height: usize,
) -> Result<GroupBBox, String> {
    // specify initial location of group bottom left
    let bottom_left = V2d {
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

fn raw_params(
    rng: &mut Pcg32,
    randomization: bool,
    rac_bounds: (Radians, Radians),
    rho_bounds: (Radians, Radians),
) -> RawParameters {
    let rac_marks = mark_between_angles(rac_bounds);
    let init_rac = RgtpDistribution::new(
        distrib_gens::specific_random(rng, 0.1, rac_marks),
        distrib_gens::random(rng, 0.1),
    );

    let rho_marks = mark_between_angles(rho_bounds);
    let init_rho = RgtpDistribution::new(
        distrib_gens::specific_random(rng, 0.1, rho_marks),
        distrib_gens::random(rng, 0.1),
    );

    defaults::RAW_PARAMS
        .modify_randomization(randomization)
        .modify_init_rac(init_rac)
        .modify_init_rho(init_rho)
}

#[allow(clippy::too_many_arguments)]
fn make_cell_group(
    rng: &mut Pcg32,
    char_quants: &CharQuantities,
    randomization: bool,
    rac_bounds: (Radians, Radians),
    rho_bounds: (Radians, Radians),
    bot_left: (Length, Length),
    num_cells: usize,
    box_width: usize,
    box_height: usize,
) -> CellGroup {
    let raw_params =
        raw_params(rng, randomization, rac_bounds, rho_bounds);
    let parameters = raw_params.refine(char_quants);
    CellGroup {
        num_cells,
        layout: group_bbox(
            num_cells,
            char_quants,
            bot_left,
            box_width,
            box_height,
        )
        .unwrap(),
        parameters,
    }
}

/// Define the cell groups that will exist in this experiment.
fn make_cell_groups(
    rng: &mut Pcg32,
    char_quants: &CharQuantities,
    randomization: bool,
    sep_in_cell_diams: usize,
) -> Vec<CellGroup> {
    let group_zero = make_cell_group(
        rng,
        char_quants,
        randomization,
        (*RAD_0, *RAD_PI),
        (*RAD_PI, *RAD_0),
        (Length(0.0), Length(0.0)),
        1,
        1,
        1,
    );
    let group_one = make_cell_group(
        rng,
        char_quants,
        randomization,
        (*RAD_PI, *RAD_0),
        (*RAD_0, *RAD_PI),
        (
            Length(0.0),
            defaults::CELL_DIAMETER.scale(sep_in_cell_diams as f64),
        ),
        1,
        1,
        1,
    );

    vec![group_zero, group_one]
}

pub fn generate(
    dirs: Directories,
    args: ExperimentArgs,
) -> Vec<Experiment> {
    let ExperimentArgs {
        toml_name,
        ty,
        final_t,
        cil_mag,
        coa_mag,
        cal_mag,
        adh_scale,
        snap_period,
        max_on_ram,
        seeds,
        int_opts,
    } = args;

    let sep_in_cell_diams =
        if let ExperimentType::Pair { sep_in_cell_diams } = &ty {
            *sep_in_cell_diams
        } else {
            panic!("Expected a Pair experiment, but got: {:?}", ty)
        };

    seeds
        .iter()
        .map(|&seed| {
            let (mut rng, randomization) = match seed {
                Some(s) => (Pcg32::seed_from_u64(s), true),
                None => (Pcg32::from_entropy(), false),
            };

            let char_quants = *defaults::CHAR_QUANTS;
            let raw_world_params = *defaults::RAW_WORLD_PARAMS;
            raw_world_params.modify_interactions(
                RawInteractionParams {
                    coa: coa_mag.map(|mag| {
                        RAW_COA_PARAMS_WITH_ZERO_MAG.modify_mag(mag)
                    }),
                    chem_attr: None,
                    bdry: None,
                    phys_contact: RawPhysicalContactParams {
                        range: RawCloseBounds {
                            zero_at: defaults::PHYS_CLOSE_DIST
                                .scale(2.0),
                            one_at: *defaults::PHYS_CLOSE_DIST,
                        },
                        adh_mag: adh_scale
                            .map(|x| defaults::ADH_MAG.scale(x)),
                        cal_mag,
                        cil_mag,
                    },
                },
            );
            let world_params = raw_world_params.refine(&char_quants);
            let cgs = make_cell_groups(
                &mut rng,
                &char_quants,
                randomization,
                sep_in_cell_diams,
            );

            let seed_string = seed
                .map_or(String::from("None"), |s| format!("{}", s));

            Experiment {
                ty: ExperimentType::NCells,
                name: format!("{}_seed={}", toml_name, seed_string),
                final_t,
                char_quants,
                world_params,
                cell_groups: cgs,
                rng,
                seed,
                snap_period,
                max_on_ram,
                int_opts,
                out_dir: (&dirs.out).clone(),
                py_main: None,
            }
        })
        .collect()
}
