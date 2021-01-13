use crate::AppState;
use druid::{
    BoxConstraints, Color, Data, Env, Event, EventCtx, LayoutCtx,
    LifeCycle, LifeCycleCtx, PaintCtx, Point, Rect, RenderContext,
    Size, UpdateCtx, Widget,
};

#[derive(Copy, Clone, Data)]
pub struct Animator {
    pub size: Size,
}

impl Animator {
    pub fn new() -> Animator {
        Animator {
            size: Size::new(800.0, 600.0),
        }
    }
}

impl Widget<AppState> for Animator {
    fn event(
        &mut self,
        _ctx: &mut EventCtx,
        _event: &Event,
        _data: &mut AppState,
        _env: &Env,
    ) {
    }

    fn lifecycle(
        &mut self,
        _ctx: &mut LifeCycleCtx,
        _event: &LifeCycle,
        _data: &AppState,
        _env: &Env,
    ) {
    }

    fn update(
        &mut self,
        _ctx: &mut UpdateCtx,
        _old_data: &AppState,
        _data: &AppState,
        _env: &Env,
    ) {
    }

    fn layout(
        &mut self,
        _ctx: &mut LayoutCtx,
        _bc: &BoxConstraints,
        _data: &AppState,
        _env: &Env,
    ) -> Size {
        self.size
    }

    fn paint(
        &mut self,
        ctx: &mut PaintCtx,
        data: &AppState,
        _env: &Env,
    ) {
        ctx.fill(
            Rect::from_origin_size(Point::new(0.0, 0.0), self.size),
            &Color::BLUE,
        )
    }
}
