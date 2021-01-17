#![allow(clippy::too_many_arguments)]
//! The entry point.
// use rust_ncc::world::hardio::Format;
use rust_ncc::{experiments, world, DEFAULT_OUTPUT_DIR};
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let exp = experiments::n_cells::generate(Some(3), 16);

    let output_dir = PathBuf::from(DEFAULT_OUTPUT_DIR);
    let mut w = world::World::new(exp, None, 10, 100);

    let now = Instant::now();
    w.simulate(3.0 * 3600.0, true);

    println!("Simulation complete. {} s.", now.elapsed().as_secs());
}
