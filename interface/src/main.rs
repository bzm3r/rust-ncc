mod animator;
mod scene;
mod view;

use crate::scene::Scene;
use crate::view::build_ui;
use druid::{AppLauncher, Data, WindowDesc};

#[derive(Clone, Data, Default)]
pub struct AppState {
    scene: Scene,
}

impl AppState {
    pub fn new() -> AppState {
        AppState {
            scene: Scene::new(),
        }
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
