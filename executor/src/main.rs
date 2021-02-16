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
    let directories = Directories::try_from(&cfg_path).unwrap();

    let exp_tomls: Vec<String> = parsed_args
        .values_of("experiments")
        .unwrap()
        .map(|s| s.into())
        .collect();

    let mut exp_toml_args = vec![];
    for exp_toml in exp_tomls.iter() {
        let fp: PathBuf = [
            &directories.exp,
            &PathBuf::from(format!("{}.toml", exp_toml)),
        ]
        .iter()
        .collect();
        exp_toml_args.push(ExperimentArgs::try_from(&fp).unwrap());
    }

    for exp_args in exp_toml_args {
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
