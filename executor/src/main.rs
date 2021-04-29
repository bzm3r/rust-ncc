// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use clap::{App, AppSettings, Arg};
use simulator::exp_setup::exp_parser::ExperimentArgs;
use simulator::{exp_setup, world, Directories};
use std::convert::TryFrom;
use std::env::current_dir;
use std::path::PathBuf;
use std::time::Instant;

pub const EXP_DIR: &str = "B:\\rust-ncc\\experiments\\";

fn main() {
    let parsed_args = App::new("simulator executor")
        .version("0.1")
        .about("Execute simulations from the command line.")
        .arg(
            Arg::with_name("config")
                .short("c")
                .long("cfg")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("experiments")
                .short("e")
                .long("exps")
                .required(true)
                .takes_value(true)
                .multiple(true)
                .min_values(1),
        )
        .setting(AppSettings::TrailingVarArg)
        .get_matches();

    let default_cfg_path: PathBuf =
        [current_dir().unwrap(), PathBuf::from("cfg")]
            .iter()
            .collect();
    let cfg_path = parsed_args
        .value_of("config")
        .map_or_else(|| default_cfg_path, PathBuf::from);
    let directories = Directories::try_from(&cfg_path).map_or_else(
        |e| panic!("{}: {:?}", &cfg_path.display(), e),
        |d| d,
    );

    let exp_jsons: Vec<String> = parsed_args
        .values_of("experiments")
        .unwrap()
        .map(|s| s.into())
        .collect();

    let mut exp_json_args = vec![];
    for exp_json in exp_jsons.iter() {
        let fp: PathBuf = [
            &directories.exp,
            &PathBuf::from(format!("{}.json", exp_json)),
        ]
        .iter()
        .collect();
        exp_json_args.push(ExperimentArgs::try_from(&fp).unwrap());
    }

    for exp_args in exp_json_args {
        let exps = exp_setup::generate(directories.clone(), exp_args);

        for exp in exps {
            let mut w = world::World::new(exp);

            let now = Instant::now();
            w.simulate(true);

            println!(
                "Simulation complete. {} s.",
                now.elapsed().as_secs()
            );
        }
    }
}
