use crate::cell::Cell;
use crate::math::v2d::V2D;
use crate::parameters::{CharQuantities, Parameters};
use crate::utils::circ_ix_plus;
use crate::world::hardio::{load_compact, Format};
use crate::world::{History, Snapshot};
use crate::NVERTS;
use lazy_static::lazy_static;
use nannou::color::gradient::Gradient;
use nannou::color::LinSrgba;
use nannou::event::WindowEvent::{KeyPressed, MousePressed};
use nannou::event::{Key, MouseButton, WindowEvent};
use nannou::geom::{Point2, Rect, Vector2};

use nannou::{color, LoopMode};
use nannou::{App, Draw, Frame};
use std::cmp::Ordering;
use std::path::PathBuf;
use std::str::FromStr;

lazy_static! {
    static ref BLACK: LinSrgba = LinSrgba::new(0.0, 0.0, 0.0, 1.0);
    static ref BLUE: LinSrgba = LinSrgba::new(0.0, 0.0, 1.0, 1.0);
    static ref RED: LinSrgba = LinSrgba::new(1.0, 0.0, 0.0, 1.0);
    static ref BLACK_TO_BLUE: Gradient<LinSrgba> =
        Gradient::new(vec![
            LinSrgba::new(0.0, 0.0, 0.0, 1.0),
            LinSrgba::new(0.0, 0.0, 1.0, 1.0)
        ]);
    static ref BLACK_TO_RED: Gradient<LinSrgba> =
        Gradient::new(vec![
            LinSrgba::new(0.0, 0.0, 0.0, 1.0),
            LinSrgba::new(1.0, 0.0, 0.0, 1.0)
        ]);
}

pub fn animate() {
    nannou::app(model).run();
}

impl From<V2D> for Vector2 {
    fn from(v: V2D) -> Self {
        Vector2 {
            x: v.x as f32,
            y: v.y as f32,
        }
    }
}

impl From<Point2<f32>> for V2D {
    fn from(v: Point2) -> Self {
        V2D { x: v.x, y: v.y }
    }
}

// fn map_rgtp_act_to_gradient(
//     rel_rgtp_act: f32,
//     rac_act: f32,
//     rho_act: f32,
//     rgtp_scale: f32,
// ) -> Gradient<LinSrgba> {
//     match rel_rgtp_act.partial_cmp(&0.0_f32).unwrap() {
//         Ordering::Less => Gradient::new(vec![
//             BLACK_TO_RED.get(rho_act / rgtp_scale),
//             LinSrgba::new(0.0, 0.0, 0.0, 0.0),
//         ]),
//         Ordering::Greater => Gradient::new(vec![
//             BLACK_TO_BLUE.get(rac_act / rgtp_scale),
//             LinSrgba::new(0.0, 0.0, 0.0, 0.0),
//         ]),
//         Ordering::Equal => Gradient::new(vec![
//             BLACK.clone(),
//             LinSrgba::new(0.0, 0.0, 0.0, 0.0),
//         ]),
//     }
// }

fn map_rgtp_act_to_color(
    rel_rgtp_act: f32,
    rac_act: f32,
    rho_act: f32,
    rgtp_scale: f32,
) -> LinSrgba {
    match rel_rgtp_act.partial_cmp(&0.0_f32).unwrap() {
        Ordering::Less => BLACK_TO_RED.get(rho_act / rgtp_scale),
        Ordering::Greater => BLACK_TO_BLUE.get(rac_act / rgtp_scale),
        Ordering::Equal => *BLACK,
    }
}

fn draw_crosshair(
    view_center: V2D,
    win_rect: Rect<f32>,
    draw: &Draw<f32>,
) {
    let weight = 2.0;
    draw.line()
        .stroke_weight(weight)
        .color(LinSrgba::new(0.0, 0.0, 0.0, 0.5))
        .points(
            Point2 {
                x: view_center.x,
                y: win_rect.bottom(),
            },
            Point2 {
                x: view_center.x,
                y: win_rect.top(),
            },
        );
    draw.line()
        .stroke_weight(weight)
        .color(LinSrgba::new(0.0, 0.0, 0.0, 0.5))
        .points(
            Point2 {
                x: win_rect.left(),
                y: view_center.y,
            },
            Point2 {
                x: win_rect.right(),
                y: view_center.y,
            },
        );
}

fn draw_cell(
    scale: f32,
    translation: V2D,
    cell: &Cell,
    rgtp_scale: f32,
    parameters: &Parameters,
    draw: &Draw<f32>,
) {
    let weight = 2.0;
    let relative_rgtp_activity =
        cell.core.calc_relative_rgtp_activity(parameters);
    let edge_vecs = cell
        .geom
        .unit_edge_vecs
        .iter()
        .zip(cell.geom.edge_lens.iter())
        .map(|(&p, &f)| f * p)
        .collect::<Vec<V2D>>();

    // let gradients = (0..NVERTS)
    //     .map(|ix| {
    //         let rel_act = relative_rgtp_activity[ix];
    //         map_rgtp_act_to_gradient(
    //             rel_act,
    //             cell.core.rac_acts[ix],
    //             cell.core.rho_acts[ix],
    //             rgtp_scale,
    //         )
    //     })
    //     .collect::<Vec<Gradient<LinSrgba>>>();

    // let colored_vertices = (0..NVERTS)
    //     .flat_map(|ix| {
    //         let v = cell.core.poly[ix];
    //         let e = edge_vecs[ix];
    //         let next_ix = circ_ix_plus(ix, NVERTS);
    //         let (g0, g1) = (&gradients[ix], &gradients[next_ix]);
    //
    //         (0..21).map(move |k| {
    //             let t = k as f32 / 20_f32;
    //             let (c0, c1) = (g0.get(t), g1.get(1.0 - t));
    //             (t * e + v, c0.mix(&c1, t))
    //         })
    //     })
    //     .collect::<Vec<(V2D, LinSrgba)>>();

    let edge_gradients = (0..NVERTS)
        .map(|ix| {
            let next_ix = circ_ix_plus(ix, NVERTS);
            Gradient::new(vec![
                map_rgtp_act_to_color(
                    relative_rgtp_activity[ix],
                    cell.core.rac_acts[ix],
                    cell.core.rho_acts[ix],
                    rgtp_scale,
                ),
                map_rgtp_act_to_color(
                    relative_rgtp_activity[next_ix],
                    cell.core.rac_acts[next_ix],
                    cell.core.rho_acts[next_ix],
                    rgtp_scale,
                ),
            ])
        })
        .collect::<Vec<Gradient<LinSrgba>>>();

    let colored_vertices = (0..NVERTS)
        .flat_map(|ix| {
            let v = cell.core.poly[ix];
            let e = edge_vecs[ix];
            let g = &edge_gradients[ix];
            (0..21).map(move |k| {
                let t = k as f32 / 20.0_f32;
                (((v - translation) + t * e).scale(scale), g.get(t))
            })
        })
        .collect::<Vec<(V2D, LinSrgba)>>();

    // Draw the polyline as a stroked path.
    draw.polyline()
        .weight(weight)
        .join_round()
        .points_colored(colored_vertices);
}

/// The application state
struct Model {
    ix: usize,
    cell_params: Vec<Parameters>,
    snapshots: Vec<Snapshot>,
    char_quants: CharQuantities,
    scale: f32,
    translation: V2D,
    crosshairs: V2D,
    draw_crosshairs: bool,
}

impl Model {
    pub fn new() -> Model {
        let out_dir = PathBuf::from_str("./output").unwrap();
        let History {
            char_quants,
            cell_params,
            snapshots,
            ..
        } = load_compact(&out_dir, Format::Bincode, "pair");
        Model {
            ix: 0,
            char_quants,
            cell_params,
            snapshots,
            scale: 1.0,
            translation: V2D::zeros(),
            crosshairs: V2D::zeros(),
            draw_crosshairs: true,
        }
    }
}

fn progress_drawing(model: &mut Model, delta: isize) {
    let new_ix = model.ix as isize + delta;
    if new_ix >= model.snapshots.len() as isize {
        model.ix = new_ix as usize - model.snapshots.len();
    } else if new_ix < 0 {
        model.ix = (model.snapshots.len() as isize + new_ix) as usize;
    } else {
        model.ix = new_ix as usize;
    }
}

fn event(app: &App, model: &mut Model, event: WindowEvent) {
    match event {
        KeyPressed(key) => key_pressed(model, key),
        MousePressed(mouse_button) => {
            mouse_pressed(app, model, mouse_button)
        }
        _ => {}
    };
}

fn key_pressed(model: &mut Model, key: Key) {
    match key {
        Key::M => progress_drawing(model, 10),
        Key::N => progress_drawing(model, -10),
        Key::X => update_view_zone(model, 0.25),
        Key::Z => update_view_zone(model, -0.25),
        Key::H => model.draw_crosshairs = !model.draw_crosshairs,
        Key::R => {
            model.crosshairs = V2D::zeros();
            model.translation = V2D::zeros();
            model.scale = 1.0;
            update_view_zone(model, 0.0);
        }
        _ => {}
    };
}

fn set_crosshairs(model: &mut Model, position: Point2<f32>) {
    model.crosshairs = position.into();
}

fn update_view_zone(model: &mut Model, delta_scale: f32) {
    model.scale = (model.scale + delta_scale).max(0.25);
    model.translation =
        model.translation + model.crosshairs.scale(1.0 / model.scale);
    model.crosshairs = V2D::zeros();
}

fn mouse_pressed(
    app: &App,
    model: &mut Model,
    mouse_button: MouseButton,
) {
    if let MouseButton::Left = mouse_button {
        set_crosshairs(model, app.mouse.position())
    };
}

/// Nannou app model
fn model(app: &App) -> Model {
    app.new_window().event(event).view(view).build().unwrap();
    app.set_loop_mode(LoopMode::Wait);

    Model::new()
}

// /// Nannou app update
// fn update(_app: &App, model: &mut Model, _update: Update) {
//     model.update();
// }

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    let win_rect = app.window_rect();
    // Draw model
    draw.background().color(color::WHITE);
    if model.draw_crosshairs {
        draw_crosshair(model.crosshairs, win_rect, &draw)
    }
    model.snapshots[model.ix]
        .cells
        .states
        .iter()
        .enumerate()
        .for_each(|(ci, c)| {
            draw_cell(
                model.scale,
                model.translation,
                c,
                model.char_quants.frac_rgtp / 10.0,
                &model.cell_params[ci],
                &draw,
            )
        });
    // Render frame
    draw.to_frame(&app, &frame).unwrap();
}
