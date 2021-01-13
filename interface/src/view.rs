use crate::AppState;
use druid::widget::Button;
use druid::{Command, FileDialogOptions, FileSpec, Target, Widget};

pub fn build_ui() -> impl Widget<AppState> {
    let binc = FileSpec::new("Bincode", &["binc"]);
    let open_dialog_options = FileDialogOptions::new()
        .allowed_types(vec![binc])
        .name_label("Simulation file")
        .title("Load simulation data")
        .force_starting_directory("../output")
        .button_text("Load");

    let open_button =
        Button::new("Open").on_click(move |ctx, _, _| {
            ctx.submit_command(Command::new(
                druid::commands::SHOW_OPEN_PANEL,
                open_dialog_options.clone(),
                Target::Auto,
            ))
        });

    open_button
    // let canvas = Canvas::default();
    //
    // let mut app_widget = Flex::column();
    // app_widget.add_child(open_button);
    // app_widget.add_child(canvas);
    // app_widget
}
