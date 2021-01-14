mod animator;
mod scene;
mod view;
mod delegate;

use crate::scene::Scene;
use crate::view::build_ui;
use druid::{AppLauncher, Data, WindowDesc};
use rust_ncc::world::History;
use std::sync::Arc;

#[derive(Clone, Data, Default)]
pub struct AppState {
    scene: Scene,
    sim_history: Arc<History>,
}

pub fn main() {
    let main_window = WindowDesc::new(build_ui)
        .title("Rust NCC")
        .window_size((800.0, 600.0));

    let app_state = AppState::default();

    AppLauncher::with_window(main_window)
        .launch(app_state)
        .expect("Failed to launch application");
}
