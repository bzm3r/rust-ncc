use crate::reader::{AsyncReader, FetchResult};
use crate::AppState;
use druid::{
    commands, AppDelegate, Command, DelegateCtx, Env, ExtEventSink,
    Handled, Selector, Target,
};

pub const FETCH_FULFILL: Selector<FetchResult> =
    Selector::new("fragment.fetch-fulfill");

pub struct Delegate {
    pub event_sink: ExtEventSink,
}

impl AppDelegate<AppState> for Delegate {
    fn command(
        &mut self,
        _ctx: &mut DelegateCtx,
        _target: Target,
        cmd: &Command,
        app: &mut AppState,
        _env: &Env,
    ) -> Handled {
        let open_command = cmd.get(commands::OPEN_FILE);
        if let Some(file_info) = open_command {
            app.reader =
                AsyncReader::spawn_for_path(file_info.path());
            return Handled::Yes;
        }
        Handled::No
    }
}
