mod animator;
mod delegate;
mod reader;
mod scene;
mod view;

use crate::delegate::Delegate;
use crate::reader::{AppContact, AsyncReader, Msg};
use crate::scene::Scene;
use crate::view::build_ui;
use crossbeam::channel::{bounded, Receiver, Sender};
use druid::im::Vector;
use druid::{AppLauncher, Data, PaintCtx, Rect, Vec2, WindowDesc};
use log::info;
use rust_ncc::hardio::WorldSnapshot;
use rust_ncc::world::py_compare::Snapshot;
use rust_ncc::world::{Snapshot, WorldInfo};
use std::sync::Arc;

pub struct ReadChannel {
    tx_to_reader: Sender<Msg>,
    rx_from_reader: Receiver<Msg>,
}

#[derive(Clone, Data)]
pub struct AppState {
    scene: Scene,
    reader_tx: Arc<Sender<AppContact>>,
    snap_offset: usize,
    world_info: Arc<WorldInfo>,
    snapshots: Arc<Vector<WorldSnapshot>>,
}

impl AppState {
    pub fn paint_scene(
        &self,
        ctx: &mut PaintCtx,
        canvas: Rect,
        zoom: f64,
        translation: Vec2,
    ) {
        self.scene.draw_snapshot(
            ctx,
            self.snapshots.get(self.snap_offset),
            &self.world_info.cell_params,
            canvas,
            zoom,
            translation,
        );
    }
}

pub fn main() {
    let main_window = WindowDesc::new(build_ui)
        .title("Rust NCC")
        .window_size((800.0, 600.0));
    let launcher = AppLauncher::with_window(main_window);
    let event_sink = launcher.get_external_handle();

    let delegate = delegate::Delegate {
        event_sink: event_sink.clone(),
    };

    let app_state = AppState {
        scene: Default::default(),
        reader_tx: ,
        snap_offset: 0,
        world_info: Arc::new(Default::default()),
        snapshots: Arc::new(Default::default()),
    };

    simple_logger::SimpleLogger::new()
        .init()
        .expect("Failed to initialize logger.");
    info!("Successfully initialized logger.");

    launcher
        .delegate(delegate)
        .launch(app_state)
        .expect("Failed to launch application");
}
