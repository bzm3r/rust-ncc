mod rkdp5;
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
mod utils;
mod world;

use std::path::PathBuf;
use crate::experiment::load_experiment;

fn main() {
    let exp = load_experiment("2020-JUL-11-test");
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\output\\");
    let mut w = world::World::new(exp, output_dir);
    w.simulate(60.0);
    w.save_history();
}
