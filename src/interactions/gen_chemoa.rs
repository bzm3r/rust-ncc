use crate::math::geometry::Poly;
use crate::math::v2d::V2D;
use crate::parameters::ChemAttrParams;
use crate::NVERTS;

#[derive(Clone)]
pub struct ChemAttrGenerator {
    center_mag: f32,
    slope: f32,
    center: V2D,
}

impl ChemAttrGenerator {
    pub fn new(params: ChemAttrParams) -> ChemAttrGenerator {
        let ChemAttrParams {
            center,
            center_mag,
            slope,
        } = params;
        ChemAttrGenerator {
            center_mag,
            slope,
            center,
        }
    }

    pub fn generate(
        &self,
        cell_polys: &[Poly],
    ) -> Vec<[f32; NVERTS]> {
        cell_polys
            .iter()
            .map(|poly| {
                let mut x_chemoas = [0.0f32; NVERTS];
                poly.verts.iter().zip(x_chemoas.iter_mut()).for_each(
                    |(&v, x)| {
                        let r = self.center_mag
                            * self.slope
                            * (v - self.center).mag();
                        if r > 0.0 {
                            *x = r
                        }
                    },
                );
                x_chemoas
            })
            .collect()
    }
}
