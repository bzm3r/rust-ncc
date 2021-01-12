use druid::{AppLauncher, WindowDesc};

mod app_state;
use app_state::AppState;

mod scene;
mod sim_data;
mod view;

use view::build_ui;

pub fn main() {
    let main_window = WindowDesc::new(build_ui)
        .title("Todo Tutorial")
        .window_size((400.0, 400.0));

    let initial_state = AppState::new();

    AppLauncher::with_window(main_window)
        .launch(initial_state)
        .expect("Failed to launch application");
}
