use crate::scene::Scene;
use druid::{
    BoxConstraints, Data, Env, Event, EventCtx, LayoutCtx, LifeCycle,
    LifeCycleCtx, PaintCtx, Size, UpdateCtx, Widget,
};

#[derive(Copy, Clone, Data)]
pub struct Animator {}

impl Animator {}

impl Widget<Scene> for Animator {
    fn event(
        &mut self,
        _ctx: &mut EventCtx,
        _event: &Event,
        _data: &mut Scene,
        _env: &Env,
    ) {
    }

    fn lifecycle(
        &mut self,
        _ctx: &mut LifeCycleCtx,
        _event: &LifeCycle,
        _data: &Scene,
        _env: &Env,
    ) {
    }

    fn update(
        &mut self,
        _ctx: &mut UpdateCtx,
        _old_data: &Scene,
        _data: &Scene,
        _env: &Env,
    ) {
    }

    fn layout(
        &mut self,
        _ctx: &mut LayoutCtx,
        _bc: &BoxConstraints,
        _data: &Scene,
        _env: &Env,
    ) -> Size {
        Size::new(800.0, 450.0)
    }

    fn paint(
        &mut self,
        _ctx: &mut PaintCtx,
        _data: &Scene,
        _env: &Env,
    ) {
    }
}
