use crate::cell::chemistry::RgtpDistribution;
use crate::exp_setup::defaults::RAW_COA_PARAMS_WITH_ZERO_MAG;
use crate::exp_setup::exp_parser::ExperimentArgs;
use crate::exp_setup::{
    defaults, CellGroup, Experiment, ExperimentType, GroupBBox,
    RgtpDistribDefs,
};
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
) -> Result<GroupBBox, String> {
    // specify initial location of group centroid
    let centroid = V2d {
        x: char_quants.normalize(&Length(0.0)),
        y: char_quants.normalize(&Length(0.0)),
    };
    let side_len = (num_cells as f64).sqrt();
    let r = GroupBBox {
        width: side_len.ceil() as usize,
        height: (num_cells as f64 / side_len).ceil() as usize,
        bottom_left: centroid,
    };
    if r.width * r.height < num_cells {
        Err(String::from(
            "Group layout area is too small to contain required number of cells.",
        ))
    } else {
        Ok(r)
    }
}

fn raw_params(
    rng: &mut Pcg32,
    rgtp_distrib_defns: &RgtpDistribDefs,
    randomization: bool,
) -> RawParameters {
    let RgtpDistribDefs { rac, rho } = rgtp_distrib_defns;

    let init_rac = RgtpDistribution::new(
        rac.acts.into_distrib(rng),
        rac.inacts.into_distrib(rng),
    );
    let init_rho = RgtpDistribution::new(
        rho.acts.into_distrib(rng),
        rho.inacts.into_distrib(rng),
    );

    defaults::RAW_PARAMS
        .modify_randomization(randomization)
        .modify_init_rac(init_rac)
        .modify_init_rho(init_rho)
}

/// Define the cell groups that will exist in this experiment.
fn make_cell_groups(
    rng: &mut Pcg32,
    char_quants: &CharQuantities,
    num_cells: usize,
    rgtp_distrib_defns: &RgtpDistribDefs,
    randomization: bool,
) -> Vec<CellGroup> {
    vec![CellGroup {
        num_cells,
        layout: group_bbox(num_cells, char_quants).unwrap(),
        parameters: raw_params(
            rng,
            rgtp_distrib_defns,
            randomization,
        )
        .refine(char_quants),
    }]
}

pub fn generate(
    dirs: Directories,
    args: ExperimentArgs,
) -> Vec<Experiment> {
    let ExperimentArgs {
        file_name: toml_name,
        ty,
        final_t,
        cil_mag,
        coa_mag,
        cal_mag,
        adh_mag: adh_scale,
        one_at,
        zero_at,
        too_close_dist,
        snap_period,
        max_on_ram,
        rgtp_distrib_defs: rgtp_distribs,
        seeds,
        int_opts,
        randomization,
    } = args;

    let num_cells = if let ExperimentType::NCells { num_cells } = &ty
    {
        *num_cells
    } else {
        panic!("Expected an n_cell experiment, but got: {:?}", ty)
    };

    seeds
        .iter()
        .map(|&seed| {
            let mut rng = Pcg32::seed_from_u64(seed);

            let char_quants = *defaults::CHAR_QUANTS;
            let raw_world_params = defaults::RAW_WORLD_PARAMS
                .modify_interactions(RawInteractionParams {
                    coa: coa_mag.map(|mag| {
                        RAW_COA_PARAMS_WITH_ZERO_MAG
                            .modify_mag(mag)
                            .modify_too_close_dist(too_close_dist)
                    }),
                    chem_attr: None,
                    bdry: None,
                    phys_contact: RawPhysicalContactParams {
                        range: RawCloseBounds { zero_at, one_at },
                        adh_mag: adh_scale
                            .map(|x| defaults::ADH_MAG.scale(x)),
                        cal_mag,
                        cil_mag,
                    },
                });
            let world_params = raw_world_params.refine(&char_quants);
            let cgs = make_cell_groups(
                &mut rng,
                &char_quants,
                num_cells,
                &rgtp_distribs,
                randomization,
            );

            Experiment {
                ty: ty.clone(),
                name: format!("{}_seed={}", toml_name, seed),
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
                run_python: false,
            }
        })
        .collect()
}
