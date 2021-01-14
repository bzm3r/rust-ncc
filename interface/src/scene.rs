use druid::{Color, Data, Point, PaintCtx, RenderContext};
use druid::kurbo::{Line};
use rust_ncc::world::Snapshot;

#[derive(Clone, Data)]
pub struct Scene {
    pub bg_color: Color,
}

impl Default for Scene {
    fn default() -> Self {
        Scene { bg_color: Color::WHITE }
    }
}

impl Scene {
    pub fn draw_snapshot(&self, ctx: &mut PaintCtx, snapshot: Option<&Snapshot>) {
        let bg_rect = ctx.size().to_rect();
        ctx.fill(bg_rect, &self.bg_color);

        if let Some(_snap) = snapshot {
            let bg_rect = ctx.size().to_rect();
            ctx.fill(bg_rect, &Color::BLACK);
        }
    }
}


