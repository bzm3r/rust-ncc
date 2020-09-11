mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;
use std::path::PathBuf;
use std::time::Instant;

pub const NVERTS: usize = 16;

fn main() {
    let exp = experiments::adh2::generate();
    let output_dir = PathBuf::from("C:\\Users\\bhmer\\Desktop\\rust-ncc\\output\\");
    let mut w = world::World::new(exp);
    let now = Instant::now();
    w.simulate(3.0 * 3600.0);
    println!("done. {} s.", now.elapsed().as_secs());
    w.save_history(&output_dir);
}
