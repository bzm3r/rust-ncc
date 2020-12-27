//! The entry point.
// #[cfg(feature = "animate")]
// mod animator;
mod cell;
mod experiments;
mod interactions;
mod math;
mod parameters;
mod utils;
mod world;

// #[cfg(feature = "animate")]
// use crate::animator::create_animation;
// use crate::math::geometry::debug_2516;
use std::path::PathBuf;
use std::time::Instant;

// /// Number of vertices per model cell.
pub const NVERTS: usize = 16;

fn main() {
    let exp = experiments::pairs::generate(Some(67));
    let output_dir = PathBuf::from(format!(
        "{}\\output",
        std::env::current_dir().unwrap().to_str().unwrap()
    ));
    let mut w = world::World::new(exp, output_dir.clone());
    let now = Instant::now();
    w.simulate(3.0 * 3600.0);
    println!("Simulation complete. {} s.", now.elapsed().as_secs());
    //#[cfg(feature = "animate")]
    //create_animation(&w.history, &output_dir.join("out.mp4"));
    w.save_history();
    //println!("test result: {}", debug_2516());
}
