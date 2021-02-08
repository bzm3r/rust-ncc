//! The entry point.
use rust_ncc::{experiments, world};
use std::time::Instant;

fn main() {
    let exp = experiments::model_compare::generate(
        Some(3),
        vec![1, 1],
        false,
        None,
        0.0,
    );
    let mut w = world::py_compare::World::new(exp);

    let now = Instant::now();
    w.simulate(0.01 * 3600.0);

    println!("Simulation complete. {} s.", now.elapsed().as_secs());
}
