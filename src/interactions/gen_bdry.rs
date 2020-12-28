#[derive(Clone)]
pub struct BdryEffectGenerator {
    shape: Vec<V2D>,
    bbox: BBox,
    skip_bb_check: bool,
    mag: f32,
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

    pub fn generate(
        &self,
        cell_polys: &[[V2D; NVERTS]],
    ) -> Vec<[f32; NVERTS]> {
        cell_polys
            .iter()
            .map(|vs| {
                let mut x_bdrys = [0.0f32; NVERTS];
                vs.iter().zip(x_bdrys.iter_mut()).for_each(
                    |(v, x)| {
                        if (self.skip_bb_check
                            && is_point_in_poly_no_bb_check(
                                v,
                                &self.shape,
                            ))
                            || is_point_in_poly(
                                v,
                                &self.bbox,
                                &self.shape,
                            )
                        {
                            *x = self.mag;
                        }
                    },
                );
                x_bdrys
            })
            .collect()
    }
}
