use rust_ncc::parameters::{CharQuantities, Parameters};
use rust_ncc::world::{History, Snapshot};

pub struct Data {
    pub snap_freq: u32,
    pub snap_ix: usize,
    pub char_quants: CharQuantities,
    pub cell_params: Vec<Parameters>,
    pub snapshots: Vec<Snapshot>,
}

impl Data {
    pub fn new(data: History) -> Data {
        let History {
            char_quants,
            cell_params,
            snapshots,
            snap_freq,
            ..
        } = data;
        Data {
            snap_freq,
            snap_ix: 0,
            char_quants,
            cell_params,
            snapshots,
        }
    }
}
