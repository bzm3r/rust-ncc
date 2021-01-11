mod data;
mod draw;
mod model;
mod ui;

use nannou::color::LinSrgba;
use nannou::conrod_core::event::DoubleClick;
use nannou::conrod_core::widget::{file_navigator, FileNavigator};
use nannou::conrod_core::{input, Color, Positionable, Widget};
use nannou::event::WindowEvent::{KeyPressed, MousePressed};
use nannou::event::{Key, MouseButton, Update, WindowEvent};
use nannou::geom::Point2;
use nannou::{App, Draw, Frame, LoopMode};
use rust_ncc::parameters::{CharQuantities, Parameters};
use rust_ncc::world::hardio::load_binc_from_path;

use crate::data::Data;
use crate::draw::{draw_cell, draw_crosshair};
use crate::model::Model;

fn main() {
    nannou::app(model).update(update).run();
}

fn increment_tstep(model: &mut Model, delta: isize) {
    if let Some(data) = model.data.as_mut() {
        let new_ix = data.snap_ix as isize + delta;
        if new_ix >= data.snapshots.len() as isize {
            data.snap_ix = new_ix as usize - data.snapshots.len();
        } else if new_ix < 0 {
            data.snap_ix =
                (data.snapshots.len() as isize + new_ix) as usize;
        } else {
            data.snap_ix = new_ix as usize;
        }
    }
}

fn event(app: &App, model: &mut Model, event: WindowEvent) {
    match event {
        KeyPressed(key) => handle_key_press(model, key),
        MousePressed(mouse_button) => {
            handle_mouse_press(app, model, mouse_button)
        }
        _ => {}
    };
}

fn update(app: &App, model: &mut Model, _update: Update) {
    // Calling `set_widgets` allows us to instantiate some widgets.
}

fn handle_key_press(model: &mut Model, key: Key) {
    match key {
        Key::M => increment_tstep(model, 10),
        Key::N => increment_tstep(model, -10),
        Key::X => update_view_zone(model, 0.25),
        Key::Z => update_view_zone(model, -0.25),
        Key::H => model.draw_crosshairs = !model.draw_crosshairs,
        Key::R => {
            model.crosshairs = Point2::zero();
            model.translation = Point2::zero();
            model.scale = 1.0;
            update_view_zone(model, 0.0);
        }
        Key::Escape => {
            model.data = None;
            model.show_file_selector = true;
        }
        Key::Grave => {
            model.show_terminal = !model.show_terminal;
        }
        _ => {}
    };
}

fn update_view_zone(model: &mut Model, delta_scale: f32) {
    model.scale = (model.scale + delta_scale).max(0.25);
    model.translation += model.crosshairs / model.scale;
}

fn handle_mouse_press(
    app: &App,
    model: &mut Model,
    mouse_button: MouseButton,
) {
    if let MouseButton::Left = mouse_button {
        model.crosshairs = app.mouse.position();
    };
    model.ui.draw_if_changed();
}

/// Nannou app model
fn model(app: &App) -> Model {
    app.new_window().event(event).view(view).build().unwrap();
    app.set_loop_mode(LoopMode::RefreshSync);
    app.set_exit_on_escape(false);

    Model::new(app)
}

fn draw_data(app: &App, model: &Model, draw: &Draw) {
    let win_rect = app.window_rect();
    // Draw data
    draw.background().color(LinSrgba::new(1.0, 1.0, 1.0, 1.0));
    if model.draw_crosshairs {
        draw_crosshair(model.crosshairs, win_rect, &draw)
    }
    if let Some(data) = model.data.as_ref() {
        data.snapshots[data.snap_ix]
            .cells
            .states
            .iter()
            .enumerate()
            .for_each(|(ci, c)| {
                draw_cell(
                    model.scale,
                    model.translation,
                    c,
                    data.char_quants.frac_rgtp / 10.0,
                    &data.cell_params[ci],
                    &draw,
                )
            });
    }
}

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    // Draw data
    if model.data.is_some() {
        draw_data(app, model, &draw);
    } else {
        draw.background().color(LinSrgba::new(0.0, 0.0, 0.0, 0.0));
    };

    // Render frame
    draw.to_frame(&app, &frame).unwrap();

    // Draw the state of the `Ui` to the frame.
    if model.show_file_selector {
        model.ui.draw_to_frame(app, &frame).unwrap();
    }
}

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
