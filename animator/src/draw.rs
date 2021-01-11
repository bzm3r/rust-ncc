use lazy_static::lazy_static;
use nannou::color::{Gradient, LinSrgba};
use nannou::geom::{Point2, Rect};
use nannou::Draw;
use rust_ncc::cell::Cell;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::Parameters;
use rust_ncc::utils::circ_ix_plus;
use rust_ncc::NVERTS;
use std::cmp::Ordering;
use std::ops::Mul;

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

fn v2d_to_point2(v: V2D) -> Point2 {
    Point2 { x: v.x, y: v.y }
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

pub fn draw_crosshair(
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

pub fn draw_cell(
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

    let edge_gradients = (0..NVERTS)
        .map(|ix| {
            let next_ix = circ_ix_plus(ix, NVERTS);
            Gradient::new(vec![
                map_rgtp_act_to_color(
                    relative_rgtp_activity[ix].to_f32(),
                    cell.core.rac_acts[ix],
                    cell.core.rho_acts[ix],
                    rgtp_scale,
                ),
                map_rgtp_act_to_color(
                    relative_rgtp_activity[next_ix].to_f32(),
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
