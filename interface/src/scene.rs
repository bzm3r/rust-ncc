use crate::sim_data::SimData;
use druid::color::BLACK;
use druid::im::Vector;
use druid::kurbo::{PathEl, PathSeg};
use druid::kurbo::Line;
use druid::widget::Painter;
use druid::{Color, Data, Env, EventCtx, LinearGradient, PaintCtx, Point, RenderContext, UnitPoint, Vec2, Key, TextLayout, FontFamily, FontWeight};
use rust_ncc::cell::Cell;
use rust_ncc::interactions::Interactions;
use rust_ncc::math::v2d::V2D;
use rust_ncc::parameters::quantity::Quantity;
use rust_ncc::parameters::{
    CharQuantities, Parameters, WorldParameters,
};
use rust_ncc::utils::circ_ix_plus;
use rust_ncc::world::{History, Snapshot};
use rust_ncc::NVERTS;
use std::cmp::Ordering;
use druid::text::{RichText, Attribute};

const SHOW_RGTP: Key<bool> = Key::new("scene.show_rgtp");
const SHOW_COA: Key<bool> = Key::new("scene.show_coa");
const SHOW_CRL: Key<bool> = Key::new("scene.show_crl");

impl Default for LinearGradient {
    fn default() -> Self {
        LinearGradient::new(
            UnitPoint::BOTTOM,
            UnitPoint::BOTTOM,
            (Color::BLACK, Color::Black),
        )
    }
}

#[derive(Clone, Data, Lens)]
pub struct Poly {
    verts: [Point; NVERTS],
    edges: [Line; NVERTS],
    rgtp_edge_colors: [LinearGradient; NVERTS],
    coa_edge_colors: [LinearGradient; NVERTS],
    crl_edge_colors: [LinearGradient; NVERTS],
    rgtp_forces: [Vec2; NVERTS],
    adh_forces: [Vec2; NVERTS],
    /// Sum of Rho GTPase mediated, cytoplasmic, and edge foces.
    sum_internal_forces: [Vec2; NVERTS],
}

impl From<&V2D> for Point {
    fn from(v: &V2D) -> Self {
        Point::new(v.x as f64, v.y as f64)
    }
}

impl From<&V2D> for Vec2 {
    fn from(v: &V2D) -> Self {
        Vec2::new(v.x as f64, v.y as f64)
    }
}

impl From<&[V2D; NVERTS]> for [Point; NVERTS] {
    fn from(v: &V2D) -> Self {
        let mut r = [Point::default(); NVERTS];
        r.iter_mut().for_each(|p| *p = Point::from(v));
        r
    }
}

impl From<&[V2D; NVERTS]> for [Vec2; NVERTS] {
    fn from(v: &V2D) -> Self {
        let mut r = [Vec2::default(); NVERTS];
        r.iter_mut().for_each(|p| *p = Vec2::from(v));
        r
    }
}

pub struct Scales {
    rgtp: f32,
    force: f32,
    crl: f32,
}

impl Poly {
    pub fn new(cell: &Cell, interactions: &Interactions, scales: Scales) -> Poly {
        let verts = cell.core.poly.into();
        let mut edges = [Line::new(Point::default(), Point::default()); NVERTS];
        for ix in 0..NVERTS {
            edges[ix] = Line::new(cell.core.poly[ix].into(), cell.core.poly[circ_ix_plus(ix, NVERTS).into());
        }

        let rgtp_edge_colors = gen_rgtp_edge_colors(cell, 0.1);
        let crl_edge_colors = gen_crl_edge_colors(interactions, 0.5);
        let coa_edge_colors = gen_coa_edge_colors(interactions, 0.5);
        let rgtp_forces: [Vec2; 16] = cell.mech.rgtp_forces.into();
        let sum_internal_forces: [Vec2; 16] =
            cell.mech.sum_forces.into();
        let adh_forces: [Vec2; 16] = interactions.x_adhs.into();

        Poly {
            verts,
            edges,
            rgtp_edge_colors,
            coa_edge_colors,
            crl_edge_colors,
            rgtp_forces,
            adh_forces,
            sum_internal_forces,
        }
    }

    pub fn paint(ctx: &mut PaintCtx, data: &Poly, env: &Env) {
        if env.get(SHOW_RGTP) {
            for (&edge, edge_color) in data.edges.iter().zip(data.rgtp_edge_colors.iter()) {
                ctx.stroke(PathSeg::Line(edge), edge_color.into(), 2.0);
            };
        }

        if env.get(SHOW_COA) {
            for (&edge, edge_color) in data.edges.iter().zip(data.coa_edge_colors.iter()) {
                ctx.stroke(PathSeg::Line(edge), edge_color.into(), 2.0);
            };
        }

        if env.get(SHOW_CRL) {
            for (&edge, edge_color) in data.edges.iter().zip(data.crl_edge_colors.iter()) {
                ctx.stroke(PathSeg::Line(edge), edge_color.into(), 2.0);
            };
        }
    }
}

#[derive(Clone, Data, Lens)]
pub struct SceneOpts {
    show_rgtp: bool,
    show_crl: bool,
    show_coa: bool,
    show_rac_forces: bool,
    show_rho_forces: bool,
    show_adh_forces: bool,
    show_tot_forces: bool,
}

#[derive(Clone, Data, Lens)]
pub struct Scene {
    tstep: u32,
    time_in_secs: f32,
    hide_crosshair: bool,
    translation: Point,
    zoom: f32,
    cross_hair: Point,
    polys: Vector<Poly>,
    opts: SceneOpts,
}

impl Scene {
    pub fn new(tstep: u32, h: &History) -> Scene {
        let snapshot = &h.snapshots[tstep as usize];
        let time_in_secs = h.char_quants.t.mul_number((tstep as f32)).number();
        let mut polys = snapshot
            .cells
            .states
            .iter()
            .zip(snapshot.cells.interactions.iter())
            .map(|(cell, interactions)| Poly::new(cell, interactions))
            .collect();

        Scene {
            tstep,
            time_in_secs,
            hide_crosshair: false,
            translation: Point::default(),
            zoom: 0.0,
            cross_hair: Point::default(),
            polys,
            opts: SceneOpts {
                show_rgtp: true,
                show_crl: false,
                show_coa: false,
                show_rac_forces: false,
                show_rho_forces: false,
                show_adh_forces: false,
                show_tot_forces: false,
            },
        }
    }

    pub fn mutate(&mut self, tstep: u32, h: &History, scene_opts: SceneOpts) {
        let snapshot = &h.snapshots[tstep as usize];
        self.tstep = tstep;
        self.time_in_secs = h.char_quants.t.mul_number((tstep as f32)).number();
        self.scene_opts = scene_opts;
        self.polys = snapshot
            .cells
            .states
            .iter()
            .zip(snapshot.cells.interactions.iter())
            .map(|(cell, interactions)| Poly::new(cell, interactions))
            .collect();
    }

    pub fn paint(ctx: &mut PaintCtx, data: &Scene, env: &Env) {
        for poly in data.polys.iter() {
            poly.paint(ctx, poly, env);
        }
        let mut time = TextLayout::default();
        time.set_text(format!("time = {} min.", (data.time_in_secs / 60.0) as u32));
        time.draw(ctx, TOP_LEFT);
    }
}

fn map_to_color(
    x: f32,
    y: f32,
    relative: Option<f32>,
    less_color: Color,
    greater_color: Color,
    scale: f32,
) -> Color {
    let relative = if let Some(rel) = relative { rel } else { x - y };
    match relative.partial_cmp(&0.0_f32).unwrap() {
        Ordering::Less => Color::BLACK + (x / scale) * less_color,
        Ordering::Greater => {
            Color::BLACK + (y / scale) * greater_color
        }
        Ordering::Equal => Color::BLACK,
    }
}

fn gen_rgtp_edge_colors(
    cell: &Cell,
    scale: f32,
) -> [LinearGradient; NVERTS] {
    let relative_rgtp_activity =
        cell.core.calc_relative_rgtp_activity(parameters);
    let mut rgtp_edge_colors = DEFAULT_LIN_GRAD;
    (0..NVERTS).for_each(|ix| {
        let next_ix = circ_ix_plus(ix, NVERTS);
        rgtp_edge_colors[ix] = LinearGradient::new(
            UnitPoint::new(0.0, 0.5),
            UnitPoint::new(1.0, 0.5),
            vec![
                map_to_color(
                    cell.core.rho_acts[ix],
                    cell.core.rac_acts[ix],
                    Some(relative_rgtp_activity[ix].to_f32()),
                    Color::RED,
                    Color::BLUE,
                    scale,
                ),
                map_to_color(
                    cell.core.rho_acts[next_ix],
                    cell.core.rac_acts[next_ix],
                    Some(relative_rgtp_activity[next_ix].to_f32()),
                    Color::RED,
                    Color::BLUE,
                    scale,
                ),
            ],
        );
    });
    rgtp_edge_colors
}

fn gen_crl_edge_colors(
    interactions: &Interactions,
    scale: f32,
) -> [LinearGradient; NVERTS] {
    let mut crl_edge_colors = DEFAULT_LIN_GRAD;
    (0..NVERTS).for_each(|ix| {
        let next_ix = circ_ix_plus(ix, NVERTS);
        crl_edge_colors[ix] = LinearGradient::new(
            UnitPoint::new(0.0, 0.5),
            UnitPoint::new(1.0, 0.5),
            vec![
                map_to_color(
                    interactions.x_cils[ix],
                    interactions.x_cals[ix],
                    None,
                    Color::YELLOW,
                    Color::AQUA,
                    scale,
                ),
                map_to_color(
                    interactions.x_cals[next_ix],
                    interactions.x_cils[next_ix],
                    None,
                    Color::YELLOW,
                    Color::AQUA,
                    scale,
                ),
            ],
        )
    });
    crl_edge_colors
}

fn gen_coa_edge_colors(
    interactions: &Interactions,
    scale: f32,
) -> [LinearGradient; NVERTS] {
    let mut coa_edge_colors = DEFAULT_LIN_GRAD;
    (0..NVERTS).for_each(|ix| {
        let next_ix = circ_ix_plus(ix, NVERTS);
        coa_edge_colors[ix] = LinearGradient::new(
            UnitPoint::new(0.0, 0.5),
            UnitPoint::new(1.0, 0.5),
            vec![
                map_to_color(
                    interactions.x_coas[ix],
                    0.0,
                    None,
                    Color::BLACK,
                    Color::FUCHSIA,
                    scale,
                ),
                map_to_color(
                    interactions.x_coas[next_ix],
                    0.0,
                    None,
                    Color::BLACK,
                    Color::FUCHSIA,
                    scale,
                ),
            ],
        )
    });
    coa_edge_colors
}
