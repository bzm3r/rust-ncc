mod cell;
mod chemistry;
mod consts;
mod experiment;
mod mat2d;
mod math;
mod mechanics;
mod parameters;
mod quantity;
mod random;
mod rkdp5;
mod schema;
mod utils;
mod world;

use crate::experiment::load_experiment;
use std::path::PathBuf;

fn main() {
    let exp = load_experiment("2020-JUL-11-test");
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\output\\");
    let mut w = world::World::new(exp);
    w.simulate(6.0 * 3600.0);
    w.save_history(&output_dir);
    w.save_geom_history(&output_dir);
    w.save_mech_history(&output_dir);
}
