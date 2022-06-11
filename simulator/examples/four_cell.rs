use simulator::exp_setup::exp_parser::ExperimentArgs;
use simulator::exp_setup::{
    generate, DistribDef, ExperimentType, RgtpDistribDef, RgtpDistribDefs,
};
use simulator::parameters::quantity::{Length, Quantity, Time};
use simulator::world::{IntegratorOpts, RkOpts};
use simulator::{world, Directories};
use std::convert::TryFrom;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;

fn main() {
    let cfg_path = PathBuf::from_str("./cfg.json").unwrap();
    let directories = Directories::try_from(&cfg_path)
        .map_or_else(|e| panic!("{}: {:?}", &cfg_path.display(), e), |d| d);
    let exp_args = ExperimentArgs {
        file_name: "example_four_cell".to_string(),
        ty: ExperimentType::NCells {
            num_cells: 4,
            chem_dist: None,
            chem_mag: None,
        },
        final_t: Time(5400.0),
        char_t: Time(1.0),
        cil_mag: 60.0,
        coa_mag: Some(24.0),
        cal_mag: Some(60.0),
        adh_scale: Some(10.0),
        adh_break: None,
        zero_at: Length(2.0).micro(),
        too_close_dist: Length(2.0).micro(),
        snap_period: Time(10.0),
        max_on_ram: 1000,
        randomization: true,
        seeds: vec![7],
        int_opts: IntegratorOpts::Rkdp5(RkOpts {
            max_iters: 10000,
            atol: 1e-3,
            rtol: 1e-3,
            init_h_scale: 0.1,
        }),
        rgtp_distrib_defs: RgtpDistribDefs {
            rac: RgtpDistribDef {
                acts: DistribDef::Random { frac: 0.1 },
                inacts: DistribDef::Random { frac: 0.3 },
            },
            rho: RgtpDistribDef {
                acts: DistribDef::Random { frac: 0.1 },
                inacts: DistribDef::Random { frac: 0.3 },
            },
        },
        crl_one_at: Default::default(),
    };
    let exp = generate(directories, exp_args)[0].clone();
    let mut w = world::World::new(exp);

    let now = Instant::now();
    w.simulate(true);

    println!("Simulation complete. {} s.", now.elapsed().as_secs());
}
