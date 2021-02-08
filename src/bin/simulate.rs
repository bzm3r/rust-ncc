//! The entry point.
use rust_ncc::{experiments, world, DEFAULT_OUTPUT_DIR};
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let seed = 3;
    println!("seed: {}", seed);
    let exp = experiments::separated_pair::generate(Some(seed));

    let mut w = world::World::new(
        exp,
        Some(PathBuf::from(DEFAULT_OUTPUT_DIR)),
        10,
        1000,
    );

    let now = Instant::now();
    w.simulate(3.0 * 3600.0, true);

    println!("Simulation complete. {} s.", now.elapsed().as_secs());
}
