use crate::AppState;
use druid::{BoxConstraints, Data, Env, Event, EventCtx, LayoutCtx, LifeCycle, LifeCycleCtx, PaintCtx, Size, UpdateCtx, Widget, KeyEvent, Code};

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
        ctx: &mut EventCtx,
        event: &Event,
        app: &mut AppState,
        _env: &Env,
    ) {
        match event {
            Event::KeyDown(KeyEvent { code, ..}) => {
                match code {
                    Code::KeyM => {
                        app.frame = {
                            let new_frame = app.frame as isize + 1;
                            if new_frame < app.sim_history.snapshots.len() as isize {
                                new_frame
                            } else {
                                0
                            }
                        } as usize;
                        ctx.request_paint();
                    },
                    Code::KeyN => {
                        app.frame = {
                            let new_frame = app.frame as isize - 1;
                            if new_frame < 0 {
                                app.sim_history.snapshots.len() - 1
                            } else {
                                new_frame as usize
                            }
                        } as usize;
                        ctx.request_paint();
                    },
                    _ => {}
                }
            },
            _ => {}
        }
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
        old_state: &AppState,
        new_state: &AppState,
        _env: &Env,
    ) {
        if !old_state.same(new_state) {
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
