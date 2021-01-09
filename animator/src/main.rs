use lazy_static::lazy_static;
use nannou::color::{Gradient, LinSrgba};
use nannou::event::WindowEvent::{KeyPressed, MousePressed};
use nannou::event::{Key, MouseButton, Update, WindowEvent};
use nannou::geom::{Point2, Rect};
use nannou::ui::prelude::*;
use nannou::{App, Draw, Frame, LoopMode};
use rust_ncc::cell::Cell;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::{CharQuantities, Parameters};
use rust_ncc::utils::circ_ix_plus;
use rust_ncc::world::hardio::{load_compact, Format};
use rust_ncc::world::{History, Snapshot};
use rust_ncc::NVERTS;
use std::cmp::Ordering;
use std::ops::Mul;
use std::path::{Path, PathBuf};
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

widget_ids! {
    struct Ids {
        out_dir,
        file_path,
    }
}

fn main() {
    nannou::app(model).run();
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

fn v2d_to_point2(v: V2D) -> Point2 {
    Point2 { x: v.x, y: v.x }
}

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
    crosshair: Point2,
    win_rect: Rect<f32>,
    draw: &Draw<f32>,
) {
    let weight = 2.0;
    draw.line()
        .stroke_weight(weight)
        .color(LinSrgba::new(0.0, 0.0, 0.0, 0.5))
        .points(
            Point2 {
                x: crosshair.x,
                y: win_rect.bottom(),
            },
            Point2 {
                x: crosshair.x,
                y: win_rect.top(),
            },
        );
    draw.line()
        .stroke_weight(weight)
        .color(LinSrgba::new(0.0, 0.0, 0.0, 0.5))
        .points(
            Point2 {
                x: win_rect.left(),
                y: crosshair.y,
            },
            Point2 {
                x: win_rect.right(),
                y: crosshair.y,
            },
        );
}

fn draw_cell(
    scale: f32,
    translation: Point2,
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
            let v = v2d_to_point2(cell.core.poly[ix]);
            let e = v2d_to_point2(edge_vecs[ix]);
            let g = &edge_gradients[ix];
            (0..21).map(move |k| {
                let t = k as f32 / 20.0_f32;
                (((v - translation) + e.mul(t)).mul(scale), g.get(t))
            })
        })
        .collect::<Vec<(Point2, LinSrgba)>>();

    // Draw the polyline as a stroked path.
    draw.polyline()
        .weight(weight)
        .join_round()
        .points_colored(colored_vertices);
}

struct Data {
    snap_ix: usize,
    char_quants: CharQuantities,
    cell_params: Vec<Parameters>,
    snapshots: Vec<Snapshot>,
}

/// The application state
struct Model {
    data: Option<Data>,
    scale: f32,
    translation: Point2,
    crosshairs: Point2,
    draw_crosshairs: bool,
    ui: Ui,
    ids: Ids,
    file_path: Option<PathBuf>,
}

impl Model {
    pub fn new(app: &App) -> Model {
        // let out_dir = PathBuf::from_str("./output").unwrap();
        // let History {
        //     char_quants,
        //     cell_params,
        //     snapshots,
        //     ..
        // } = load_compact(&out_dir, Format::Bincode, "pair");
        // Create the UI.
        let mut ui = app.new_ui().build().unwrap();

        // Generate some ids for our widgets.
        let ids = Ids::new(ui.widget_id_generator());

        Model {
            scale: 1.0,
            translation: Point2::zero(),
            crosshairs: Point2::zero(),
            draw_crosshairs: true,
            ui,
            ids,
            file_path: None,
            data: None,
        }
    }
}

fn progress_drawing(model: &mut Model, delta: isize) {
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
        KeyPressed(key) => key_pressed(model, key),
        MousePressed(mouse_button) => {
            mouse_pressed(app, model, mouse_button)
        }
        _ => {}
    };
}

fn update(_app: &App, model: &mut Model, _update: Update) {
    // Calling `set_widgets` allows us to instantiate some widgets.
    let ui = &mut model.ui.set_widgets();

    fn file_navigator<'a>() -> widget::FileNavigator<'a> {
        widget::FileNavigator::directories(Path::new("../"))
    }
}

fn key_pressed(model: &mut Model, key: Key) {
    match key {
        Key::M => progress_drawing(model, 10),
        Key::N => progress_drawing(model, -10),
        Key::X => update_view_zone(model, 0.25),
        Key::Z => update_view_zone(model, -0.25),
        Key::H => model.draw_crosshairs = !model.draw_crosshairs,
        Key::R => {
            model.crosshairs = Point2::zero();
            model.translation = Point2::zero();
            model.scale = 1.0;
            update_view_zone(model, 0.0);
        }
        _ => {}
    };
}

fn update_view_zone(model: &mut Model, delta_scale: f32) {
    model.scale = (model.scale + delta_scale).max(0.25);
    model.translation += model.crosshairs / model.scale;
}

fn mouse_pressed(
    app: &App,
    model: &mut Model,
    mouse_button: MouseButton,
) {
    if let MouseButton::Left = mouse_button {
        model.crosshairs = app.mouse.position();
    };
}

/// Nannou app model
fn model(app: &App) -> Model {
    app.new_window().event(event).view(view).build().unwrap();
    app.set_loop_mode(LoopMode::Wait);

    Model::new(app)
}

// /// Nannou app update
// fn update(_app: &App, model: &mut Model, _update: Update) {
//     model.update();
// }

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
        Some(model.crosshairs)
    } else {
        None
    };

    // Render frame
    draw.to_frame(&app, &frame).unwrap();
}
