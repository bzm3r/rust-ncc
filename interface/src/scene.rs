use druid::{Color, Data};

#[derive(Clone, Data)]
pub struct Scene {
    pub color: Color,
}

impl Default for Scene {
    fn default() -> Self {
        Scene { color: Color::WHITE }
    }
}
