mod cell;
mod experiment;
mod math;
mod parameters;
mod utils;
mod world;

use crate::experiment::load_experiment;
use std::path::PathBuf;

pub const NVERTS: usize = 16;

fn main() {
    let exp = load_experiment("2020-JUL-11-test");
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\output\\");
    let mut w = world::World::new(exp);
    w.simulate(3.0 * 3600.0);
    w.save_history(&output_dir);
    w.save_geom_history(&output_dir);
    w.save_mech_history(&output_dir);
}
