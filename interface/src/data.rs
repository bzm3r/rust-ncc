use druid::{Data, Env, EventCtx};
use rust_ncc::parameters::quantity::Quantity;
use rust_ncc::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use rust_ncc::world::Snapshot;

#[derive(Clone, Data, Lens)]
pub struct AppState {
    terminal: Terminal,
    simulation: Option<Simulation>,
    drawing: Drawing,
}

impl AppState {
    pub fn load(_ctx: &mut EventCtx, data: &mut Self, _env: &Env) {
        unimplemented!();
    }

    pub fn click_term_in(
        _ctx: &mut EventCtx,
        data: &mut Self,
        _env: &Env,
    ) {
        unimplemented!()
    }

    pub fn input_term_in(
        _ctx: &mut EventCtx,
        data: &mut Self,
        _env: &Env,
    ) {
        unimplemented!()
    }
}

#[derive(Clone, Data, Lens)]
pub struct Terminal {
    u_in: UserInput,
    out: TermOut,
}

pub type UserInput = String;
pub type TermOut = String;

#[derive(Clone, PartialEq, Lens)]
pub struct Simulation {
    tstep: u32,
    snap_freq: u32,
    pub char_quants: CharQuantities,
    pub world_params: WorldParameters,
    pub cell_params: Vec<Parameters>,
    //TODO: make snapshot (ultimately `Cells`) use `im` vec?
    // what about serialization/deserialization into CBOR?
    snapshots: Vec<Snapshot>,
}

impl Data for Simulation {
    fn same(&self, other: &Self) -> bool {
        self.tstep == other.tstep
            && self.snap_freq == other.snap_freq
            && self.char_quants == other.char_quants
            && self.world_params == other.world_params
            && self
                .cell_params
                .iter()
                .zip(other.cell_params.iter())
                .all(|(s, o)| s == 0)
            && self
                .snapshots
                .iter()
                .zip(other.snapshots.iter())
                .all(|(s, o)| s == o)
    }
}

pub type Seconds = f32;

impl Simulation {
    pub fn state_around(&self, time: Seconds) -> Snapshot {
        let t = (time / self.char_quants.t.number()).floor() as u32;
        self.state_at(t / self.snap_freq)
    }

    pub fn state_at(&self, tstep: u32) -> Snapshot {
        self.snapshot[tstep]
    }
}

#[derive(Clone, Data, Lens)]
pub struct Drawing {
    hide_crosshair: bool,
    translation: [f32; 2],
    zoom: f32,
    cross_hair: [f32; 2],
    // should cell poly related stuff go here?
}
