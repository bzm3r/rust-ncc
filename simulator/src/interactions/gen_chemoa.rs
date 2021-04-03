use crate::math::geometry::Poly;
use crate::math::min_f64s;
use crate::math::v2d::V2d;
use crate::parameters::ChemAttrParams;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, Serialize)]
pub struct ChemAttrGenerator {
    center_mag: f64,
    slope: f64,
    center: V2d,
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

    pub fn update(&self) {}

    pub fn generate(
        &self,
        cell_polys: &[Poly],
    ) -> Vec<[f64; NVERTS]> {
        cell_polys
            .iter()
            .map(|poly| {
                let mut x_chemoas = [0.0f64; NVERTS];
                poly.verts.iter().zip(x_chemoas.iter_mut()).for_each(
                    |(&v, x_chemoa)| {
                        let chemo_signal = self.center_mag
                            * (1.0
                                - self.slope
                                    * (v - self.center).mag());
                        if chemo_signal > 0.0 {
                            *x_chemoa = chemo_signal
                        }
                    },
                );
                let min_chemoa = min_f64s(&x_chemoas);
                (0..NVERTS).for_each(|ix| {
                    x_chemoas[ix] -= min_chemoa;
                });
                x_chemoas
            })
            .collect()
    }
}
