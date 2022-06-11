// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use clap::{AppSettings, Arg, Command};
use simulator::exp_setup::exp_parser::ExperimentArgs;
use simulator::{exp_setup, world, Directories};
use std::convert::TryFrom;
use std::path::PathBuf;
use std::time::Instant;

pub const EXP_DIR: &str = "B:\\rust-ncc\\experiments\\";

fn main() {
    let parsed_args = Command::new("simulator executor")
        .version("0.1")
        .about("Execute simulations from the command line.")
        .arg(
            Arg::new("config")
                .short('c')
                .long("cfg")
                .required(true)
                .help("Configuration file. A default config can be found in ./cfg.json")
                .takes_value(true),
        )
        .arg(
            Arg::new("experiments")
                .short('e')
                .long("exps")
                .required(true)
                .takes_value(true)
                .multiple_occurrences(true)
                .help("Experiment(s) to be run. Example experiments can be found in ./example-experiments")
                .min_values(1),
        )
        .get_matches();

    let cfg_path = parsed_args.value_of("config").map(PathBuf::from).unwrap();

    let directories = Directories::try_from(&cfg_path);
    if let Err(e) = directories {
        panic!(
            "Error loading configuration file ({}): {:?}",
            &cfg_path.display(),
            e
        );
    }
    let directories = directories.unwrap();

    let exp_jsons: Vec<PathBuf> = parsed_args
        .values_of("experiments")
        .unwrap()
        .map(Into::<String>::into)
        .map(PathBuf::from)
        .collect();

    let mut exp_json_args = vec![];
    for fp in exp_jsons.iter() {
        let args = ExperimentArgs::try_from(fp);
        match args {
            Ok(args) => exp_json_args.push(args),
            Err(e) => {
                panic!("Failed to load experiment ({}): {:?}", fp.display(), e)
            }
        }
    }

    for exp_args in exp_json_args {
        let exps = exp_setup::generate(directories.clone(), exp_args);

        for exp in exps {
            let mut w = world::World::new(exp);

            let now = Instant::now();
            w.simulate(true);

            println!("Simulation complete. {} s.", now.elapsed().as_secs());
        }
    }
}
