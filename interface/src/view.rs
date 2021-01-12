use druid::{widget::Label, Widget};

use crate::app_state::*;

pub fn build_ui() -> impl Widget<AppState> {
    Label::new("Hello")
}
