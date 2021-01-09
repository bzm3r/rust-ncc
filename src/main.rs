#![allow(clippy::too_many_arguments)]
//! The entry point.
mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;

use crate::world::hardio::Format;
use std::path::PathBuf;
use std::time::Instant;

/// Number of vertices per model cell.
pub const DEFAULT_OUTPUT_DIR: &str = "./output";

fn main() {
    let exp = experiments::pair::generate(Some(3));

    let output_dir = PathBuf::from(DEFAULT_OUTPUT_DIR);
    let mut w = world::World::new(exp, output_dir.clone());

    let now = Instant::now();
    w.simulate(3.0 * 3600.0, 5);
    println!("Simulation complete. {} s.", now.elapsed().as_secs());
    let now = Instant::now();
    w.save_history(true, vec![Format::Bincode]).unwrap();
    println!(
        "Finished saving history. {} s.",
        now.elapsed().as_secs()
    );
}
