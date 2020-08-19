mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;
use std::path::PathBuf;

pub const NVERTS: usize = 16;

fn main() {
    let exp = experiments::single::generate();
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\output\\");
    let mut w = world::World::new(exp);
    w.simulate(3.0 * 3600.0);
    w.save_history(&output_dir);
}
