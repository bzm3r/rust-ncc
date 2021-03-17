use crate::animator::Animator;
use crate::AppState;
use druid::widget::{Button, Flex};
use druid::{Command, FileDialogOptions, FileSpec, Target, Widget};

pub fn build_ui() -> impl Widget<AppState> {
    let binc = FileSpec::new("Bincode", &["binc"]);
    let open_dialog_options = FileDialogOptions::new()
        .allowed_types(vec![binc])
        .name_label("Simulation file")
        .title("Load simulation data")
        .button_text("Load");

    let open_button =
        Button::new("Open").on_click(move |ctx, _, _| {
            ctx.submit_command(Command::new(
                druid::commands::SHOW_OPEN_PANEL,
                open_dialog_options.clone(),
                Target::Auto,
            ))
        });

    let animator = Animator::new();
    Flex::column().with_child(open_button).with_flex_child(animator, 1.0)
}
