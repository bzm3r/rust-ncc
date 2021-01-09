use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::dat_sym4d::SymCcVvDat;
use crate::interactions::generate_contacts;
use crate::math::geometry::{BBox, LineSeg2D, Poly};
use crate::parameters::CoaParams;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct VertexPairInfo {
    dist: f32,
    num_intersects: f32,
}

impl VertexPairInfo {
    pub fn infinity() -> VertexPairInfo {
        VertexPairInfo {
            dist: f32::INFINITY,
            num_intersects: f32::INFINITY,
        }
    }
}

#[derive(Clone, Deserialize, Serialize)]
pub struct CoaGenerator {
    dat: SymCcVvDat<VertexPairInfo>,
    contact_bbs: Vec<BBox>,
    contacts: SymCcDat<bool>,
    params: CoaParams,
}

// pub fn self_intersects(
//     vi: u32,
//     poly: &[V2D; NVERTS],
//     lseg: &LineSeg2D,
// ) -> bool {
//     let ui = circ_ix_minus(vi, NVERTS);
//     for i in 0..NVERTS {
//         if !(i == vi || i == ui)
//             && lseg
//                 .intersects_lseg(&LineSeg2D::new(
//                     &poly[i],
//                     &poly[circ_ix_plus(i, NVERTS)],
//                 ))
//                 .is_some()
//         {
//             return true;
//         }
//     }
//     false
// }

// pub struct VertexInfo {
//     ci: u32,
//     v: V2D,
// }

/// Calculate clearance and distance.
pub fn calc_pair_info(
    ci: usize,
    oci: usize,
    lseg: LineSeg2D,
    poly_a: &Poly,
    poly_b: &Poly,
    cell_polys: &[Poly],
) -> VertexPairInfo {
    if !(lseg.intersects_poly(poly_a) || lseg.intersects_poly(poly_b))
    {
        let dist = lseg.mag();
        let clearance = cell_polys
            .iter()
            .enumerate()
            .map(|(pi, poly)| {
                if pi != ci && pi != oci && lseg.intersects_poly(poly)
                {
                    1.0
                } else {
                    0.0
                }
            })
            .sum::<f32>();
        return VertexPairInfo {
            dist,
            num_intersects: clearance,
        };
    }
    VertexPairInfo::infinity()
}

impl CoaGenerator {
    /// Calculates a matrix storing whether two vertices have clear line of sight if in contact range.
    pub fn new(
        cell_polys: &[Poly],
        params: CoaParams,
    ) -> CoaGenerator {
        let num_cells = cell_polys.len();
        let contact_bbs = cell_polys
            .iter()
            .map(|cp| cp.bbox.expand_by(params.range))
            .collect::<Vec<BBox>>();
        let contacts = generate_contacts(&contact_bbs);
        let mut dat =
            SymCcVvDat::empty(num_cells, VertexPairInfo::infinity());

        for (ci, poly) in cell_polys.iter().enumerate() {
            for (ocj, opoly) in
                cell_polys[(ci + 1)..].iter().enumerate()
            {
                let oci = (ci + 1) + ocj;
                for (vi, v) in poly.verts.iter().enumerate() {
                    for (ovi, ov) in opoly.verts.iter().enumerate() {
                        dat.set(
                            ci,
                            vi,
                            oci,
                            ovi,
                            calc_pair_info(
                                ci,
                                oci,
                                LineSeg2D::new(v, ov),
                                poly,
                                opoly,
                                cell_polys,
                            ),
                        )
                    }
                }
            }
        }
        CoaGenerator {
            dat,
            contact_bbs,
            contacts,
            params,
        }
    }

    pub fn update(&mut self, ci: usize, cell_polys: &[Poly]) {
        let this_poly = cell_polys[ci];
        let bb = this_poly.bbox.expand_by(self.params.range);
        self.contact_bbs[ci] = bb;
        // Update contacts.
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contacts.set(ci, oci, obb.intersects(&bb))
            }
        }
        for (oci, other_poly) in cell_polys.iter().enumerate() {
            if oci == ci || !self.contacts.get(ci, oci) {
                continue;
            }
            for (vi, v) in this_poly.verts.iter().enumerate() {
                for (ovi, ov) in other_poly.verts.iter().enumerate() {
                    self.dat.set(
                        ci,
                        vi,
                        oci,
                        ovi,
                        calc_pair_info(
                            ci,
                            oci,
                            LineSeg2D::new(v, ov),
                            &this_poly,
                            other_poly,
                            cell_polys,
                        ),
                    );
                }
            }
        }
    }

    pub fn generate(&self) -> Vec<[f32; NVERTS]> {
        let num_cells = self.contacts.num_cells;
        let mut all_x_coas = vec![[0.0f32; NVERTS]; num_cells];
        let CoaParams {
            los_penalty,
            mag,
            distrib_exp,
            ..
        } = self.params;
        for (ci, x_coas) in all_x_coas.iter_mut().enumerate() {
            for (vi, x_coa) in x_coas.iter_mut().enumerate() {
                for oci in 0..num_cells {
                    if oci != ci {
                        for ovi in 0..NVERTS {
                            let VertexPairInfo {
                                dist,
                                num_intersects,
                            } = self.dat.get(ci, vi, oci, ovi);
                            *x_coa = mag * (distrib_exp * dist).exp()
                                / (num_intersects + 1.0)
                                    .powf(los_penalty);
                        }
                    }
                }
            }
        }
        all_x_coas
    }
}
