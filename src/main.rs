mod experiments;
mod interactions;
mod math;
mod model_cell;
mod parameters;
mod utils;
mod world;
use std::path::PathBuf;
use std::time::Instant;

pub const NVERTS: usize = 16;

fn main() {
    let exp = experiments::single::generate();
    let output_dir = PathBuf::from(format!(
        "{}\\output",
        std::env::current_dir().unwrap().to_str().unwrap()
    ));
    let mut w = world::World::new(exp);
    let now = Instant::now();
    w.simulate(3.0 * 3600.0);
    println!("done. {} s.", now.elapsed().as_secs());
    w.save_history(&output_dir);
}
