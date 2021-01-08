use crate::cell::Cell;
use crate::math::v2d::V2D;
use crate::parameters::{Parameters, WorldParameters};
use crate::utils::circ_ix_plus;
use crate::world::hardio::{load_compact, Format};
use crate::world::{MiniHistory, MiniSnapshot};
use crate::NVERTS;
use nannou::color;
use nannou::color::gradient::Gradient;
use nannou::color::Rgb;
use nannou::event::Key;
use nannou::geom::Vector2;
use nannou::{App, Draw, Frame};
use std::path::PathBuf;
use std::str::FromStr;
use tracing::field::display;

pub fn animate() {
    nannou::app(model).run();
}

impl From<&V2D> for Vector2 {
    fn from(v: &V2D) -> Self {
        Vector2 {
            x: v.x as f32,
            y: v.y as f32,
        }
    }
}

fn display_cell(
    cell: &Cell,
    parameters: &Parameters,
    draw: &Draw<f32>,
) {
    let weight = 2.0;
    let relative_rgtp_activity =
        cell.core.calc_relative_rgtp_activity(parameters);
    let colors = (0..NVERTS)
        .map(|vi| {
            match (
                relative_rgtp_activity[vi] > 0.0,
                relative_rgtp_activity[circ_ix_plus(vi, NVERTS)]
                    > 0.0,
            ) {
                (true, true) => {
                    Gradient::new(&[Rgb::new(0.0, 0.0, 1.0)])
                }
                (true, false) => Gradient::new(&[
                    Rgb::new(0.0, 0.0, 1.0),
                    Rgb::new(1.0, 0.0, 0.0),
                ]),
                (false, true) => Gradient::new(&[
                    Rgb::new(1.0, 0.0, 0.0),
                    Rgb::new(0.0, 0.0, 1.0),
                ]),
                (false, false) => {
                    Gradient::new(&[Rgb::new(1.0, 0.0, 0.0)])
                }
            }
        })
        .collect::<Vec<Gradient<Rgb>>>();
    let colored_vertices = cell
        .core
        .poly
        .iter()
        .zip(colors.iter())
        .collect::<Vec<(&V2D, &Gradient<Rgb>)>>();

    // Draw the polyline as a stroked path.
    draw.polyline()
        .weight(weight)
        .join_round()
        .points_colored_closed(colored_vertices);
}

/// The application state
struct Model {
    ix: usize,
    world_params: WorldParameters,
    cell_params: Vec<Parameters>,
    snapshots: Vec<MiniSnapshot>,
}

impl Model {
    pub fn new() -> Model {
        let out_dir = PathBuf::from_str("./output").unwrap();
        let MiniHistory {
            world_params,
            cell_params,
            snapshots,
        } = load_compact(&out_dir, Format::Bincode, "pair");
        Model {
            ix: 0,
            world_params,
            cell_params,
            snapshots,
        }
    }
}

fn key_pressed(_app: &App, model: &mut Model, key: Key) {
    let delta = match key {
        Key::M => 10,
        Key::N => -10,
        _ => 0,
    };
    let new_ix = model.ix as isize + delta;
    if new_ix > model.data.len() as isize {
        model.ix = new_ix as usize - model.data.len();
    } else if new_ix < 0 {
        model.ix = (model.data.len() as isize + new_ix) as usize;
    } else {
        model.ix = new_ix as usize;
    }
}

/// Nannou app model
fn model(app: &App) -> Model {
    app.new_window()
        .key_pressed(key_pressed)
        .view(view)
        .build()
        .unwrap();

    Model::new()
}

// /// Nannou app update
// fn update(_app: &App, model: &mut Model, _update: Update) {
//     model.update();
// }

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    // Draw model
    draw.background().color(color::WHITE);
    model.snapshots[model.ix]
        .cells
        .states
        .iter()
        .enumerate()
        .for_each(|(ci, c)| {
            display_cell(c, &model.cell_params[ci], &draw)
        });
    // Render frame
    draw.to_frame(&app, &frame).unwrap();
}
