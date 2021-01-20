mod animator;
mod delegate;
mod reader;
mod scene;
mod view;

use crate::delegate::Delegate;
use crate::reader::{AsyncReader, Msg};
use crate::scene::Scene;
use crate::view::build_ui;
use crossbeam::channel::{bounded, Receiver, Sender};
use druid::im::Vector;
use druid::{AppLauncher, Data, PaintCtx, Rect, Vec2, WindowDesc};
use log::info;
use rust_ncc::world::{Snapshot, WorldInfo};
use std::sync::Arc;
use std::thread;

pub struct ReadChannel {
    tx_to_reader: Sender<Msg>,
    rx_from_reader: Receiver<Msg>,
}
#[derive(Clone, Data)]
pub struct AppState {
    scene: Scene,
    read_channel: Arc<ReadChannel>,
    snap_offset: usize,
    world_info: Arc<WorldInfo>,
    snapshots: Arc<Vector<Snapshot>>,
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

    let (tx_to_app, rx_from_reader): (Sender<Msg>, Receiver<Msg>) =
        bounded(1);
    let (tx_to_reader, reader) = AsyncReader::new(tx_to_app);
    thread::spawn(move || {
        let mut reader = reader;
        reader.work_loop();
    });

    let app_state = AppState {
        scene: Default::default(),
        read_channel: Arc::new(ReadChannel {
            tx_to_reader,
            rx_from_reader,
        }),
        snap_offset: 0,
        world_info: Arc::new(Default::default()),
        snapshots: Arc::new(Default::default()),
    };

    simple_logger::SimpleLogger::new()
        .init()
        .expect("Failed to initialize logger.");
    info!("Successfully initialized logger.");

    AppLauncher::with_window(main_window)
        .delegate(Delegate)
        .launch(app_state)
        .expect("Failed to launch application");
}
