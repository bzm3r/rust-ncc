use crate::math::geometry::{calc_poly_area, BBox};
use crate::math::v2d::V2D;
use crate::NVERTS;

pub struct CellPoly {
    bbox: BBox,
    verts: [V2D; NVERTS],
}

impl CellPoly {
    pub fn area(&self) -> f32 {
        calc_poly_area(&self.verts)
    }

    pub fn expand_bbox(&self, range: f32) -> BBox {
        self.bbox.expand_by(range)
    }
}
