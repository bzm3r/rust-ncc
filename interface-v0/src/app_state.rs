use crate::canvas::Canvas;
use druid::{Data, Lens};
use rust_ncc::world::History;
use std::sync::Arc;

#[derive(Clone, Data, Lens)]
pub struct AppState {
    //terminal: Terminal,
    canvas: Arc<Canvas>,
    pub sim_data: Arc<History>,
}

impl AppState {
    pub fn new() -> AppState {
        AppState {
            //terminal: Terminal {},
            canvas: Arc::new(Canvas::default()),
            sim_data: Arc::new(History::default()),
        }
    }
}

// #[derive(Clone, Data, Lens)]
// pub struct Terminal {
//     u_in: UserInput,
//     out: TermOut,
// }
//
// pub type UserInput = String;
// pub type TermOut = String;



