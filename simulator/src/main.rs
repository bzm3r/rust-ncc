/*
 * // Copyright © 2022 Brian Merchant.
 * //
 * // Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
 * // http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
 * // <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
 * // option. This file may not be copied, modified, or distributed
 * // except according to those terms.
 */

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

use tracing::{info, Level};

fn main() {
    tracing_subscriber::fmt()
        .compact()
        // all spans/events with a level higher than TRACE (e.g, info, warn, etc.)
        // will be written to stdout.
        .with_max_level(Level::TRACE)
        // sets this to be the default, global collector for this application.
        .init();

    let cfg_path = PathBuf::from_str("B:/rust-ncc/cfg.json").unwrap();
    let directories = Directories::try_from(&cfg_path)
        .map_or_else(|e| panic!("{}: {:?}", &cfg_path.display(), e), |d| d);
    let exp_args = ExperimentArgs {
        file_name: "perftest_16_cell".to_string(),
        ty: ExperimentType::NCells {
            num_cells: 16,
            chem_dist: None,
            chem_mag: None,
        },
        final_t: Time(540.0),
        char_t: Time(0.1),
        cil_mag: 60.0,
        coa_mag: Some(14.0),
        cal_mag: Some(60.0),
        adh_scale: None,
        adh_break: None,
        zero_at: Length(2.0).micro(),
        too_close_dist: Length(2.0).micro(),
        snap_period: Time(1.0),
        max_on_ram: 1000,
        randomization: true,
        seeds: vec![7],
        int_opts: IntegratorOpts::Rkdp5(RkOpts {
            max_iters: 1000,
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

    info!("Simulation complete. {} s.", now.elapsed().as_secs());
}
