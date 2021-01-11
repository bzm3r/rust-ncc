use crate::data::Data;
use crate::model::Model;
use nannou::conrod_core::event::DoubleClick;
use nannou::conrod_core::widget::{file_navigator, FileNavigator};
use nannou::conrod_core::{input, Color, Positionable, Widget};
use rust_ncc::world::hardio::load_binc_from_path;

pub fn init_ui(model: &mut Model) {
    if model.show_file_selector {
        let mut ui = model.ui.set_widgets();
        for event in
            FileNavigator::with_extension(&model.directory, &["binc"])
                .top_left_with_margin(20.0)
                .set(model.ids.file_selector, &mut ui)
        {
            if let file_navigator::Event::DoubleClick(
                DoubleClick {
                    button: input::MouseButton::Left,
                    ..
                },
                paths,
            ) = event
            {
                model.file_path = paths[0].clone();
                model.data = Some(Data::new(load_binc_from_path(
                    &model.file_path,
                )));
                model.show_file_selector = false;
            }
        }
    } else {
        model.ui.clear_with(Color::Rgba(0.0, 0.0, 0.0, 0.0))
    }
}
