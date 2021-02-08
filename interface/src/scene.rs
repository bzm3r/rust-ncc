use druid::kurbo::{Line, PathSeg};
use druid::{
    Affine, Color, Data, LinearGradient, PaintCtx, Point, Rect,
    RenderContext, UnitPoint, Vec2,
};
use rust_ncc::cell::Cell;
use rust_ncc::interactions::RelativeRgtpActivity;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::Parameters;
use rust_ncc::utils::circ_ix_plus;
use rust_ncc::world::Snapshot;
use rust_ncc::NVERTS;

#[derive(Clone, Data)]
pub struct Scene {
    pub bg_color: Color,
}

impl Default for Scene {
    fn default() -> Self {
        Scene {
            bg_color: Color::WHITE,
        }
    }
}

pub fn line_from_v2ds(u: V2D, v: V2D) -> Line {
    Line::new(
        Point {
            x: u.x as f64,
            y: u.y as f64,
        },
        Point {
            x: v.x as f64,
            y: v.y as f64,
        },
    )
}

impl Scene {
    pub fn draw_snapshot(
        &self,
        ctx: &mut PaintCtx,
        snapshot: Option<&Snapshot>,
        _parameters: &[Parameters],
        canvas: Rect,
        zoom: f64,
        translation: Vec2,
    ) {
        ctx.fill(canvas, &self.bg_color);

        if let Some(snap) = snapshot {
            let coordinate_transform = Affine::FLIP_Y
                * Affine::translate(Vec2::new(
                    0.5 * canvas.width(),
                    -0.5 * canvas.height(),
                ))
                * Affine::scale(zoom)
                * Affine::translate(translation);
            let translation =
                coordinate_transform * Affine::translate(translation);
            ctx.transform(Affine::scale(zoom) * translation);
            ctx.save().unwrap();
            for cell in snap.cells.states.iter() {
                let edge_colors =
                    colors_from_rgtp_activity(cell, 0.1);
                for (ix, color) in edge_colors.iter().enumerate() {
                    let u = cell.core.poly[ix];
                    let v = cell.core.poly[circ_ix_plus(ix, NVERTS)];
                    ctx.stroke(
                        PathSeg::Line(line_from_v2ds(u, v)),
                        color,
                        2.0,
                    );
                }
            }
            ctx.restore().unwrap();
        }
    }
}

pub fn mix_colors(w0: f64, c0: Color, w1: f64, c1: Color) -> Color {
    let (r0, g0, b0, a0) = c0.as_rgba();
    let (r1, g1, b1, a1) = c1.as_rgba();
    Color::rgba(
        w0 * r0 + w1 * r1,
        w0 * g0 + w1 * g1,
        w0 * b0 + w1 * b1,
        w0 * a0 + w1 * a1,
    )
}

// pub fn map_to_color(
//     left: f32,
//     right: f32,
//     left_color: Color,
//     right_color: Color,
//     scale: f32,
// ) -> Color {
//     match (left - right).partial_cmp(&0.0_f32) {
//         Some(order) => match order {
//             Ordering::Less => mix_colors(
//                 1.0,
//                 Color::BLACK,
//                 (left / scale) as f64,
//                 left_color,
//             ),
//             Ordering::Greater => mix_colors(
//                 1.0,
//                 Color::BLACK,
//                 (right / scale) as f64,
//                 right_color,
//             ),
//             Ordering::Equal => Color::BLACK,
//         },
//         None => panic!("map_color: could not compare f32s successfully using partial_cmp!")
//     }
// }

const MID_0: UnitPoint = UnitPoint::new(0.0, 0.5);
const MID_1: UnitPoint = UnitPoint::new(1.0, 0.5);

fn colors_from_relative_rgtp_activity(
    cell: &Cell,
    scale: f32,
    parameters: &Parameters,
) -> Vec<LinearGradient> {
    let relative_rgtp_activity =
        cell.core.calc_relative_rgtp_activity(parameters);
    let vertex_colors = relative_rgtp_activity
        .iter()
        .map(|&rel_act| match rel_act {
            RelativeRgtpActivity::RhoDominant(rho) => mix_colors(
                (rho / scale) as f64,
                Color::RED,
                1.0 - (rho / scale) as f64,
                Color::BLACK,
            ),
            RelativeRgtpActivity::RacDominant(rac) => mix_colors(
                (rac / scale) as f64,
                Color::BLUE,
                1.0 - (rac / scale) as f64,
                Color::BLACK,
            ),
        })
        .collect::<Vec<Color>>();
    let mut rgtp_edge_colors: Vec<LinearGradient> =
        Vec::with_capacity(NVERTS);
    for ix in 0..NVERTS {
        let this_color = vertex_colors[ix].clone();
        let next_ix = circ_ix_plus(ix, NVERTS);
        let next_color = vertex_colors[next_ix].clone();

        rgtp_edge_colors.push(LinearGradient::new(
            MID_0.clone(),
            MID_1.clone(),
            (this_color, next_color),
        ));
    }
    rgtp_edge_colors
}

fn colors_from_rgtp_activity(
    cell: &Cell,
    scale: f32,
) -> Vec<LinearGradient> {
    let rac_acts = cell.core.rac_acts;
    let rho_acts = cell.core.rho_acts;
    let vertex_colors = rac_acts
        .iter()
        .zip(rho_acts.iter())
        .map(|(rac_act, rho_act)| {
            let (rac_act, rho_act) =
                (rac_act / scale, rho_act / scale);
            let sum_act = rac_act + rho_act;
            let rac = rac_act / sum_act;
            let rho = rho_act / sum_act;

            mix_colors(
                rac as f64,
                Color::BLUE,
                rho as f64,
                Color::RED,
            )
        })
        .collect::<Vec<Color>>();
    let mut rgtp_edge_colors: Vec<LinearGradient> =
        Vec::with_capacity(NVERTS);
    for ix in 0..NVERTS {
        let this_color = vertex_colors[ix].clone();
        let next_ix = circ_ix_plus(ix, NVERTS);
        let next_color = vertex_colors[next_ix].clone();

        rgtp_edge_colors.push(LinearGradient::new(
            MID_0.clone(),
            MID_1.clone(),
            (this_color, next_color),
        ));
    }
    rgtp_edge_colors
}
