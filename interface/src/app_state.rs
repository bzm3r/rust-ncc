use crate::scene::Scene;
use crate::sim_data::SimData;
use druid::{Data, Env, EventCtx};
use rust_ncc::parameters::quantity::Quantity;
use rust_ncc::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use rust_ncc::world::Snapshot;

#[derive(Clone, Data, Lens)]
pub struct AppState {
    terminal: Terminal,
    picture: Scene,
}

impl AppState {
    pub fn set_pic(&mut self, sim_dat: SimData) {
        self.picture = Scene::from(sim_dat);
    }
}

#[derive(Clone, Data, Lens)]
pub struct Terminal {
    u_in: UserInput,
    out: TermOut,
}

pub type UserInput = String;
pub type TermOut = String;
