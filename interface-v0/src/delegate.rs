use crate::app_state::AppState;
use crate::artist::Artist;
use druid::{
    commands, AppDelegate, Command, DelegateCtx, Env, Handled, Target,
};
use rust_ncc::world::hardio::load_binc_from_path;
use std::sync::Arc;

pub struct Delegate;

impl AppDelegate<AppState> for Delegate {
    fn command(
        &mut self,
        _ctx: &mut DelegateCtx,
        _target: Target,
        cmd: &Command,
        data: &mut AppState,
        _env: &Env,
    ) -> Handled {
        if let Some(file_info) = cmd.get(commands::OPEN_FILE) {
            data.sim_data =
                Arc::new(load_binc_from_path(file_info.path()));
            return Handled::Yes;
        }
        Handled::No
    }
}
