mod view;
mod scene;
mod animator;

use crate::view::build_ui;
use druid::{AppLauncher, Data, WindowDesc};

#[derive(Clone, Copy, Data, Default)]
pub struct AppState {}

impl AppState {
    pub fn new() -> AppState {
        AppState {}
    }
}

pub fn main() {
    let main_window = WindowDesc::new(build_ui)
        .title("Rust NCC")
        .window_size((800.0, 600.0));

    let app_state = AppState::new();

    AppLauncher::with_window(main_window)
        .launch(app_state)
        .expect("Failed to launch application");
}
