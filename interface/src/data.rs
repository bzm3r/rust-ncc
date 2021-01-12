use druid::Data;
use rust_ncc::parameters::quantity::Quantity;
use rust_ncc::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use rust_ncc::world::Snapshot;

#[derive(Clone, Data)]
pub struct AppState {
    terminal: Terminal,
    simulation: Simulation,
    drawing: Drawing,
}

#[derive(Clone, Data)]
pub struct Terminal {
    u_in: UserInput,
    out: TermOut,
}

pub type UserInput = String;
pub type TermOut = String;

#[derive(Clone, Data)]
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

impl Data for CharQuantities {
    fn same(&self, other: &Self) -> bool {
        unimplemented!()
    }
}

impl Data for WorldParameters {
    fn same(&self, other: &Self) -> bool {
        unimplemented!()
    }
}

impl Data for Parameters {
    fn same(&self, other: &Self) -> bool {
        unimplemented!()
    }
}

impl Data for Snapshot {
    fn same(&self, other: &Self) -> bool {
        unimplemented!()
    }
}

pub type Seconds = f32;

impl Simulation {
    pub fn state_around(&self, time: Seconds) -> Snapshot {
        let t = (time / self.char_quants.t.number()).floor() as u32;
        self.state_at(t / self.snap_freq);
    }

    pub fn state_at(tstep: u32) -> Snapshot {
        self.snapshot[tstep]
    }
}

#[derive(Clone, Data)]
pub struct Drawing {}
