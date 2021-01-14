use crate::AppState;
use druid::{BoxConstraints, Color, Data, Env, Event, EventCtx, LayoutCtx, LifeCycle, LifeCycleCtx, PaintCtx, RenderContext, Size, UpdateCtx, Widget};

#[derive(Copy, Clone, Data)]
pub struct Animator {
}

impl Animator {
    pub fn new() -> Animator {
        Animator {}
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
        bc: &BoxConstraints,
        _data: &AppState,
        _env: &Env,
    ) -> Size {
        bc.max()
    }

    fn paint(
        &mut self,
        ctx: &mut PaintCtx,
        _data: &AppState,
        _env: &Env,
    ) {
        let rect = ctx.size().to_rect();
        ctx.fill(
            rect,
            &Color::BLUE,
        )
    }
}
