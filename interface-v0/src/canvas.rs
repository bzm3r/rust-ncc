// const SHOW_RGTP: Key<bool> = Key::new("scene.show_rgtp");
// const SHOW_COA: Key<bool> = Key::new("scene.show_coa");
// const SHOW_CRL: Key<bool> = Key::new("scene.show_crl");

use crate::artist::Artist;
use druid::{
    BoxConstraints, Env, Event, EventCtx, LayoutCtx, LifeCycle,
    LifeCycleCtx, PaintCtx, Point, Size, UpdateCtx, Widget,
};

#[derive(Default)]
pub struct Canvas {
    tstep: u32,
    hide_crosshair: bool,
    cross_hair: Point,
    translation: Point,
    zoom: f32,
    opts: Opts,
    artist: Artist,
}

impl Widget<Artist> for Canvas {
    fn event(
        &mut self,
        _ctx: &mut EventCtx,
        _event: &Event,
        _data: &mut Artist,
        _env: &Env,
    ) {
    }

    fn lifecycle(
        &mut self,
        _ctx: &mut LifeCycleCtx,
        _event: &LifeCycle,
        _data: &Artist,
        _env: &Env,
    ) {
    }

    fn update(
        &mut self,
        _ctx: &mut UpdateCtx,
        _old_data: &Artist,
        _data: &Artist,
        _env: &Env,
    ) {
    }

    fn layout(
        &mut self,
        _ctx: &mut LayoutCtx,
        _bc: &BoxConstraints,
        _data: &Artist,
        _env: &Env,
    ) -> Size {
        Size::new(600.0, 600.0)
    }

    fn paint(
        &mut self,
        _ctx: &mut PaintCtx,
        _data: &Artist,
        _env: &Env,
    ) {
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq)]
pub struct Opts {
    show_rgtp: bool,
    show_crl: bool,
    show_coa: bool,
    show_rac_forces: bool,
    show_rho_forces: bool,
    show_adh_forces: bool,
    show_tot_forces: bool,
}
