use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::dat_sym4d::SymCcVvDat;
use crate::interactions::generate_contacts;
use crate::math::geometry::{BBox, LineSeg2D};
use crate::math::v2d::V2D;
use crate::parameters::CoaParams;
use crate::utils::{circ_ix_minus, circ_ix_plus};
use crate::NVERTS;

#[derive(Clone, Copy)]
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

#[derive(Clone)]
pub struct CoaGenerator {
    dat: SymCcVvDat<VertexPairInfo>,
    contact_bbs: Vec<BBox>,
    contacts: SymCcDat<bool>,
    params: CoaParams,
}

pub fn self_intersects(
    vi: usize,
    poly: &[V2D; NVERTS],
    lseg: &LineSeg2D,
) -> bool {
    let ui = circ_ix_minus(vi, NVERTS);
    for i in 0..NVERTS {
        if !(i == vi || i == ui)
            && lseg
                .intersects_lseg(&LineSeg2D::new(
                    &poly[i],
                    &poly[circ_ix_plus(i, NVERTS)],
                ))
                .is_some()
        {
            return true;
        }
    }
    false
}

pub struct VertexInfo {
    ci: usize,
    vi: usize,
    v: V2D,
}

/// Calculate clearance and distance.
pub fn calc_pair_info(
    vert: VertexInfo,
    poly: &[V2D; NVERTS],
    overt: VertexInfo,
    opoly: &[V2D; NVERTS],
    cell_polys: &[[V2D; NVERTS]],
) -> VertexPairInfo {
    let VertexInfo { ci, vi, v } = vert;
    let VertexInfo {
        ci: oci,
        vi: ovi,
        v: ov,
    } = overt;
    let lseg = LineSeg2D::new(&v, &ov);
    if !(self_intersects(vi, poly, &lseg)
        || self_intersects(ovi, opoly, &lseg))
    {
        let dist = (v - ov).mag();
        let clearance = cell_polys
            .iter()
            .enumerate()
            .map(|(pi, poly)| {
                if pi != ci
                    && pi != oci
                    && lseg.intersects_poly(None, poly)
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
        cell_bbs: &[BBox],
        cell_polys: &[[V2D; NVERTS]],
        params: CoaParams,
    ) -> CoaGenerator {
        let num_cells = cell_bbs.len();
        let contact_bbs = cell_bbs
            .iter()
            .map(|bb| bb.expand_by(params.range))
            .collect::<Vec<BBox>>();
        let contacts = generate_contacts(&contact_bbs);
        let mut dat =
            SymCcVvDat::empty(num_cells, VertexPairInfo::infinity());

        for (ci, poly) in cell_polys.iter().enumerate() {
            for (ocj, opoly) in
                cell_polys[(ci + 1)..].iter().enumerate()
            {
                let oci = (ci + 1) + ocj;
                for (vi, v) in poly.iter().enumerate() {
                    for (ovi, ov) in opoly.iter().enumerate() {
                        dat.set(
                            ci,
                            vi,
                            oci,
                            ovi,
                            calc_pair_info(
                                VertexInfo { ci, vi, v: *v },
                                poly,
                                VertexInfo {
                                    ci: oci,
                                    vi: ovi,
                                    v: *ov,
                                },
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

    pub fn update(
        &mut self,
        ci: usize,
        bb: &BBox,
        poly: &[V2D; NVERTS],
        cell_polys: &[[V2D; NVERTS]],
    ) {
        let bb = bb.expand_by(self.params.range);
        self.contact_bbs[ci] = bb;
        // Update contacts.
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contacts.set(ci, oci, obb.intersects(&bb))
            }
        }
        for (oci, opoly) in cell_polys.iter().enumerate() {
            if oci == ci || !self.contacts.get(ci, oci) {
                continue;
            }
            for (vi, v) in poly.iter().enumerate() {
                for (ovi, ov) in opoly.iter().enumerate() {
                    self.dat.set(
                        ci,
                        vi,
                        oci,
                        ovi,
                        calc_pair_info(
                            VertexInfo { ci, vi, v: *v },
                            poly,
                            VertexInfo {
                                ci: oci,
                                vi: ovi,
                                v: *ov,
                            },
                            opoly,
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
