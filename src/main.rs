use crate::experiment::load_experiment;
mod cell;
mod mat2d;
mod math;
mod parameters;
mod quantity;
mod random;
mod utils;
mod world;
mod consts;
mod experiment;
mod chemistry;
mod mechanics;

fn main() {
    let exp = load_experiment("2020-JUL-11-test");
    let mut w = world::World::new(exp);
    w.simulate(60.0);
}
