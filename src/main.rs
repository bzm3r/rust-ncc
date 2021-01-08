#![allow(clippy::too_many_arguments)]
//! The entry point.
#[cfg(feature = "animate")]
mod animator;
mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;

#[cfg(feature = "animate")]
use crate::animator::create_animation;
use crate::world::hardio::Format;
use std::path::PathBuf;
use std::time::Instant;

// /// Number of vertices per model cell.
pub const NVERTS: usize = 16;

fn main() {
    let exp = experiments::four_cells::generate(Some(3));

    let output_dir = PathBuf::from("./output");
    let mut w = world::World::new(exp, output_dir.clone());

    let now = Instant::now();
    w.simulate(3.0 * 3600.0, 10);
    println!("Simulation complete. {} s.", now.elapsed().as_secs());
    let now = Instant::now();
    w.save_history(true, vec![Format::Cbor, Format::Bincode])
        .unwrap();
    println!(
        "Finished saving history. {} s.",
        now.elapsed().as_secs()
    );

    #[cfg(feature = "animate")]
    create_animation(&w.history, &output_dir.join("out.mp4"));
}
