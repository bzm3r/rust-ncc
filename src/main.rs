#![allow(clippy::too_many_arguments)]
//! The entry point.
//#[cfg(feature = "animate")]
mod animator;
mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;

#[cfg(feature = "animate")]
use crate::animator::animate;
use crate::world::hardio::Format;
use std::path::PathBuf;
use std::time::Instant;

/// Number of vertices per model cell.
pub const NVERTS: usize = 16;
pub const OUTPUT_DIR: &'static str = "./output";

fn main() {
    let exp = experiments::pair::generate(Some(3));

    let output_dir = PathBuf::from(OUTPUT_DIR);
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

    #[cfg(feature = "animate")]
    animate()
}
