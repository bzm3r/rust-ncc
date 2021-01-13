use crate::artist::Scales;
use crate::utils::{point_arr_from, point_from};
use druid::kurbo::{Line, PathSeg};
use druid::{
    Color, Env, LinearGradient, PaintCtx, Point, RenderContext,
    UnitPoint,
};
use rust_ncc::cell::Cell;
use rust_ncc::interactions::{Interactions, RelativeRgtpActivity};
use rust_ncc::parameters::Parameters;
use rust_ncc::utils::circ_ix_plus;
use rust_ncc::NVERTS;
use std::cmp::Ordering;

#[derive(Clone)]
pub struct Polygon {
    verts: [Point; NVERTS],
    edges: [Line; NVERTS],
    rgtp_edge_colors: Vec<LinearGradient>,
    // coa_edge_colors: Vec<LinearGradient>,
    // crl_edge_colors: Vec<LinearGradient>,
    // rgtp_forces: [Vec2; NVERTS],
    // adh_forces: [Vec2; NVERTS],
    // /// Sum of Rho GTPase mediated, cytoplasmic, and edge foces.
    // sum_internal_forces: [Vec2; NVERTS],
}

impl Polygon {
    pub fn new(
        cell: &Cell,
        interactions: &Interactions,
        parameters: &Parameters,
        scales: Scales,
    ) -> Polygon {
        let verts = point_arr_from(&cell.core.poly);
        let mut edges =
            [Line::new(Point::default(), Point::default()); NVERTS];
        for ix in 0..NVERTS {
            edges[ix] = Line::new(
                point_from(&cell.core.poly[ix]),
                point_from(&cell.core.poly[circ_ix_plus(ix, NVERTS)]),
            );
        }

        let rgtp_edge_colors = relative_rgtp_activity_coloration(
            cell,
            scales.rgtp,
            parameters,
        );
        // let crl_edge_colors = [default_lin_grad(); NVERTS];
        // let coa_edge_colors = [default_lin_grad(); NVERTS];
        // let rgtp_forces: [Vec2; 16] =
        //     vec2_arr_from(&cell.mech.rgtp_forces);
        // let sum_internal_forces: [Vec2; 16] =
        //     vec2_arr_from(&cell.mech.sum_forces);
        // let adh_forces: [Vec2; 16] =
        //     vec2_arr_from(&interactions.x_adhs);

        Polygon {
            verts,
            edges,
            rgtp_edge_colors,
            // coa_edge_colors,
            // crl_edge_colors,
            // rgtp_forces,
            // adh_forces,
            // sum_internal_forces,
        }
    }

    pub fn paint(ctx: &mut PaintCtx, data: &Polygon, env: &Env) {
        for (&edge, edge_color) in
            data.edges.iter().zip(data.rgtp_edge_colors.iter())
        {
            ctx.stroke(PathSeg::Line(edge), edge_color.into(), 2.0);
        }

        // if env.get(SHOW_COA) {
        //     for (&edge, edge_color) in
        //         data.edges.iter().zip(data.coa_edge_colors.iter())
        //     {
        //         ctx.stroke(
        //             PathSeg::Line(edge),
        //             edge_color.into(),
        //             2.0,
        //         );
        //     }
        // }
        //
        // if env.get(SHOW_CRL) {
        //     for (&edge, edge_color) in
        //         data.edges.iter().zip(data.crl_edge_colors.iter())
        //     {
        //         ctx.stroke(
        //             PathSeg::Line(edge),
        //             edge_color.into(),
        //             2.0,
        //         );
        //     }
        // }
    }
}

pub fn mix_colors(w0: f64, c0: Color, w1: f64, c1: Color) -> Color {
    let (r0, b0, g0, a0) = c0.as_rgba();
    let (r1, b1, g1, a1) = c1.as_rgba();
    Color::rgba(
        w0 * r0 + w1 * r1,
        w0 * g0 + w1 * g1,
        w0 * b0 + w1 * b1,
        w0 * a0 + w1 * a1,
    )
}

pub fn map_to_color(
    left: f32,
    right: f32,
    left_color: Color,
    right_color: Color,
    scale: f32,
) -> Color {
    match (left - right).partial_cmp(&0.0_f32) {
        Some(order) => match order {
            Ordering::Less => mix_colors(
                1.0,
                Color::BLACK,
                (left / scale) as f64,
                left_color,
            ),
            Ordering::Greater => mix_colors(
                1.0,
                Color::BLACK,
                (right / scale) as f64,
                right_color,
            ),
            Ordering::Equal => Color::BLACK,
        },
        None => panic!("map_color: could not compare f32s successfully using partial_cmp!")
    }
}

const MID_0: UnitPoint = UnitPoint::new(0.0, 0.5);
const MID_1: UnitPoint = UnitPoint::new(1.0, 0.5);

fn relative_rgtp_activity_coloration(
    cell: &Cell,
    scale: f32,
    parameters: &Parameters,
) -> Vec<LinearGradient> {
    let relative_rgtp_activity =
        cell.core.calc_relative_rgtp_activity(parameters);
    let vertex_colors =
        relative_rgtp_activity.iter().map(|&rel_act| match rel_act {
            RelativeRgtpActivity::RhoDominant(rho) => mix_colors(
                (rho / scale as f64) as f64,
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
        });
    let mut rgtp_edge_colors: Vec<LinearGradient> =
        Vec::with_capacity(NVERTS);
    for ix in 0..NVERTS {
        let this_color = vertex_colors[ix];
        let next_ix = circ_ix_plus(ix, NVERTS);
        let next_color = vertex_colors[next_ix];

        rgtp_edge_colors.push(LinearGradient::new(
            MID_0.clone(),
            MID_1.clone(),
            (this_color, next_color),
        ));
    }
    rgtp_edge_colors
}
