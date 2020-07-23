use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use std::os::raw::c_void;

#[repr(C)]
#[derive(Clone, Copy, Default)]
struct P2D {
    x: f32,
    y: f32,
}

#[repr(C)]
struct State {
    positions: [P2D; 6],
    dat: [f32; 6],
}

impl State {
    fn new() -> State {
        let distr: Uniform<f32> = Uniform::new_inclusive(0.0, 1.0);
        let mut rng = thread_rng();

        let mut positions = [P2D::default(); 6];
        let mut dat = [0.0; 6];

        for i in 0..6 {
            positions[i] = P2D {
                x: rng.sample(distr),
                y: rng.sample(distr),
            };

            dat[i] = rng.sample(distr);
        }

        State {
            positions,
            dat,
        }
    }
}

fn main() {
    let state = Box::new(State::new());
    let state_ptr = Box::into_raw(state);
    
}
