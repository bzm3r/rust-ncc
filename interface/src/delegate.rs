use crate::reader::HistoryReader;
use crate::AppState;
use druid::{
    commands, AppDelegate, Command, DelegateCtx, Env, Handled, Target,
};

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
        let open_command = cmd.get(commands::OPEN_FILE);
        if let Some(file_info) = open_command {
            app.read_channel.send()
            app.read_channel = HistoryReader::from(file_info.path());
            return Handled::Yes;
        }
        Handled::No
    }
}
