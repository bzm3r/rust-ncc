use druid::{Color, Data};

#[derive(Clone, Data)]
pub struct Scene {
    color: Color,
}

impl Scene {
    pub fn new() -> Scene {
        Scene::default()
    }
}

impl Default for Scene {
    fn default() -> Self {
        Scene { color: Color::BLUE }
    }
}
