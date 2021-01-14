use druid::{Color, Data, Point, PaintCtx, RenderContext, Rect, Affine, Vec2};
use druid::kurbo::{Line, PathSeg};
use rust_ncc::world::Snapshot;
use rust_ncc::math::v2d::V2D;
use rust_ncc::NVERTS;
use rust_ncc::utils::circ_ix_plus;

#[derive(Clone, Data)]
pub struct Scene {
    pub bg_color: Color,
}

impl Default for Scene {
    fn default() -> Self {
        Scene { bg_color: Color::WHITE }
    }
}

pub fn line_from_v2ds(u: V2D, v: V2D) -> Line {
    Line::new(Point {
        x: u.x as f64,
        y: u.y as f64,
    }, Point {
        x: v.x as f64,
        y: v.y as f64,
    })
}


impl Scene {
    pub fn draw_snapshot(&self, ctx: &mut PaintCtx, snapshot: Option<&Snapshot>, canvas: Rect) {
        ctx.fill(canvas, &self.bg_color);

        if let Some(snap) = snapshot {
            ctx.transform(Affine::FLIP_Y * Affine::translate(Vec2::new(0.5 * canvas.width(), -0.5 * canvas.height())));
            ctx.save().unwrap();
            for state in snap.cells.states.iter() {
                for ix in 0..NVERTS {
                    let u = state.core.poly[ix];
                    let v = state.core.poly[circ_ix_plus(ix, NVERTS)];
                    ctx.stroke(PathSeg::Line(line_from_v2ds(u, v)), &Color::BLACK, 2.0);
                }
            }
            ctx.restore().unwrap();
        }
    }
}


