use crate::math::geometry::{is_point_in_poly, BBox, Poly};
use crate::math::v2d::V2d;
use crate::parameters::BdryParams;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, Serialize)]
pub struct BdryEffectGenerator {
    shape: Vec<V2d>,
    bbox: BBox,
    skip_bb_check: bool,
    mag: f64,
}

impl BdryEffectGenerator {
    pub fn new(params: BdryParams) -> BdryEffectGenerator {
        let BdryParams {
            shape,
            bbox,
            skip_bb_check,
            mag,
        } = params;
        BdryEffectGenerator {
            shape,
            bbox,
            skip_bb_check,
            mag,
        }
    }

    pub fn generate(&self, cell_polys: &[Poly]) -> Vec<[f64; NVERTS]> {
        cell_polys
            .iter()
            .map(|poly| {
                let mut x_bdrys = [0.0f64; NVERTS];
                poly.verts
                    .iter()
                    .zip(x_bdrys.iter_mut())
                    .for_each(|(v, x)| {
                        let in_bdry = if self.skip_bb_check {
                            is_point_in_poly(v, None, &self.shape)
                        } else {
                            is_point_in_poly(v, Some(&self.bbox), &self.shape)
                        };
                        if in_bdry {
                            *x = self.mag;
                        }
                    });
                x_bdrys
            })
            .collect()
    }
}
