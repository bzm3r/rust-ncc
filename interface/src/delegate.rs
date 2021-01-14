use crate::AppState;
use druid::{AppDelegate, DelegateCtx, Target, Command, Env, Handled, commands};
use std::sync::Arc;
use rust_ncc::world::hardio::load_binc_from_path;

pub struct Delegate;

impl AppDelegate<AppState> for Delegate {
    fn command(
        &mut self,
        _ctx: &mut DelegateCtx,
        _target: Target,
        cmd: &Command,
        app: &mut AppState,
        _env: &Env,
    ) -> Handled {
        if let Some(file_info) = cmd.get(commands::OPEN_FILE) {
            app.sim_history =
                Arc::new(load_binc_from_path(file_info.path()));
            return Handled::Yes;
        }
        Handled::No
    }
}
