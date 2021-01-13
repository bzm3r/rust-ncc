use druid::{AppLauncher, WindowDesc};

mod app_state;
use app_state::AppState;

mod artist;
mod canvas;
mod delegate;
mod poly;
mod polygon;
mod utils;
mod view;

use delegate::Delegate;
use view::build_ui;

pub fn main() {
    let main_window = WindowDesc::new(build_ui)
        .title("Rust NCC")
        .window_size((800.0, 600.0));

    let initial_state = AppState::new();

    AppLauncher::with_window(main_window)
        .delegate(Delegate {})
        .launch(initial_state)
        .expect("Failed to launch application");
}
