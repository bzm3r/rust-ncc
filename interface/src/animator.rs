use crate::reader::{FetchRequest, Msg};
use crate::AppState;
use druid::kurbo::Line;
use druid::{
    BoxConstraints, Code, Color, Data, Env, Event, EventCtx,
    KeyEvent, LayoutCtx, LifeCycle, LifeCycleCtx, PaintCtx, Point,
    Rect, RenderContext, Size, UpdateCtx, Vec2, Widget,
};
use log::info;
use rust_ncc::world::Snapshot;

#[derive(Copy, Clone, Default, Data)]
pub struct Crosshairs {
    pos: Point,
    hidden: bool,
}

#[derive(Copy, Clone, Data)]
pub struct Animator {
    crosshairs: Crosshairs,
    zoom: f64,
    translation: Vec2,
}

impl Animator {
    pub fn new() -> Animator {
        Animator {
            crosshairs: Crosshairs::default(),
            zoom: 1.0,
            translation: Vec2::new(0.0, 0.0),
        }
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
            Event::MouseDown(mouse) => {
                if !ctx.is_focused() {
                    ctx.request_focus();
                    ctx.set_active(true);
                }

                self.translation = mouse.pos - self.crosshairs.pos;
                self.crosshairs.pos = mouse.pos;
                ctx.request_paint();
            }

            Event::KeyDown(KeyEvent { code, .. }) => {
                match code {
                    Code::KeyM => {
                        // increment scene by one frame
                        if app.snap_offset == app.snapshots.len() {
                            app.snapshots = match app
                                .read_tx
                                .tx_to_reader
                                .send(Msg::FromApp(
                                    FetchRequest::FetchNext,
                                )) {
                                Ok() => {}
                                Err(_) => {}
                            };
                            app.snap_offset = 1;
                        }
                        ctx.request_paint();
                    }
                    Code::KeyH => {
                        // toggle visibility of crosshairs
                        self.crosshairs.hidden =
                            !self.crosshairs.hidden;
                        ctx.request_paint();
                    }
                    Code::Equal => {
                        self.zoom += 0.25;
                        ctx.request_paint();
                    }
                    Code::Minus => {
                        self.zoom -= 0.25;
                        if self.zoom < 0.0 {
                            self.zoom = 0.0;
                        }
                        ctx.request_paint();
                    }
                    _ => {}
                }
            }
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
        app.paint_scene(ctx, canvas, self.zoom, self.translation);
        self.draw_crosshairs(ctx, canvas);
    }
}

impl Animator {
    fn draw_crosshairs(&self, ctx: &mut PaintCtx, canvas: Rect) {
        if !(self.crosshairs.hidden) {
            let horizontal = Line::new(
                Point::new(canvas.x0, self.crosshairs.pos.y),
                Point::new(canvas.x1, self.crosshairs.pos.y),
            );
            let vertical = Line::new(
                Point::new(
                    canvas.x0 + self.crosshairs.pos.x,
                    canvas.y0,
                ),
                Point::new(
                    canvas.x0 + self.crosshairs.pos.x,
                    canvas.y1,
                ),
            );
            ctx.stroke(
                horizontal,
                &Color::grey(0.25).with_alpha(0.3),
                1.0,
            );
            ctx.stroke(
                vertical,
                &Color::grey(0.25).with_alpha(0.3),
                1.0,
            );
        }
    }
}
