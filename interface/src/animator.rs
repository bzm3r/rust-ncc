use crate::AppState;
use druid::{BoxConstraints, Data, Env, Event, EventCtx, LayoutCtx, LifeCycle, LifeCycleCtx, PaintCtx, Size, UpdateCtx, Widget};

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
        ctx: &mut UpdateCtx,
        old_data: &AppState,
        data: &AppState,
        _env: &Env,
    ) {
        if !old_data.same(data) {
            ctx.request_paint();
        }
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
        app: &AppState,
        _env: &Env,
    ) {
        let canvas = ctx.size().to_rect();
        app.paint_scene(ctx, canvas);
    }
}
