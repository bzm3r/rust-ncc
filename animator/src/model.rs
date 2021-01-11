use crate::data::Data;
use nannou::geom::Point2;
use nannou::ui::widget_ids;
use nannou::{App, Ui};
use std::path::PathBuf;

widget_ids! {
    pub struct Ids {
        file_selector,
    }
}

/// The application state
pub struct Model {
    pub data: Option<Data>,
    pub scale: f32,
    pub translation: Point2,
    pub crosshairs: Point2,
    pub draw_crosshairs: bool,
    pub ui: Ui,
    pub ids: Ids,
    pub directory: PathBuf,
    pub file_path: PathBuf,
    pub show_file_selector: bool,
    pub show_terminal: bool,
}

impl Model {
    pub fn new(app: &App) -> Model {
        let mut ui = app.new_ui().build().unwrap();

        let ids = Ids::new(ui.widget_id_generator());

        Model {
            scale: 1.0,
            translation: Point2::zero(),
            crosshairs: Point2::zero(),
            draw_crosshairs: true,
            ui,
            ids,
            directory: PathBuf::from("../output"),
            file_path: PathBuf::from(""),
            data: None,
            show_file_selector: true,
            show_terminal: false,
        }
    }
}
