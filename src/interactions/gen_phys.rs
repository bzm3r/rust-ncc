use crate::interactions::dat4d::CvCvDat;
use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::{generate_contacts, RgtpState};
use crate::math::geometry::{seg_to_point_dist, BBox};
use crate::math::v2d::V2D;
use crate::parameters::PhysicalContactParams;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use std::cmp::Ordering;

pub type ClosestDistLoc = f32;
pub type ClosestDist = f32;

/// Generates CIL/CAL related interaction information. In other words,
/// interactions that require cells to engage in physical contact.
#[derive(Clone)]
pub struct PhysicalContactGenerator {
    ///
    dat: CvCvDat<(ClosestDistLoc, ClosestDist)>,
    pub contact_bbs: Vec<BBox>,
    pub contacts: SymCcDat<bool>,
    params: PhysicalContactParams,
}

pub struct PhysContactFactors {
    pub adh: Vec<[V2D; NVERTS]>,
    pub cil: Vec<[f32; NVERTS]>,
    pub cal: Vec<[f32; NVERTS]>,
}

impl PhysicalContactGenerator {
    /// Calculate distances between vertices of cells in contact.
    pub fn new(
        cell_bbs: &[BBox],
        cell_polys: &[[V2D; NVERTS]],
        params: PhysicalContactParams,
    ) -> PhysicalContactGenerator {
        let num_cells = cell_bbs.len();
        let mut dat =
            CvCvDat::empty(num_cells, (f32::INFINITY, f32::INFINITY));
        let contact_bbs = cell_bbs
            .iter()
            .map(|bb| bb.expand_by(params.range))
            .collect::<Vec<BBox>>();
        let contacts = generate_contacts(&contact_bbs);
        for (ci, vs) in cell_polys.iter().enumerate() {
            for (oci, ovs) in cell_polys.iter().enumerate() {
                if ci != oci && contacts.get(ci, oci) {
                    for (vi, v) in vs.iter().enumerate() {
                        for (ovi, ov) in ovs.iter().enumerate() {
                            let ov2 = &ovs[circ_ix_plus(ovi, NVERTS)];
                            let (t, d) =
                                seg_to_point_dist(v, ov, ov2);
                            match d.partial_cmp(&params.range) {
                                Some(ord) => match ord {
                                    Ordering::Greater => {
                                        continue;
                                    }
                                    _ => {
                                        dat.set(
                                            ci,
                                            vi,
                                            oci,
                                            ovi,
                                            (t, d),
                                        );
                                        break;
                                    }
                                },
                                _ => panic!(
                                    "cannot compare {} and  {}",
                                    d, &params.range
                                ),
                            }
                        }
                    }
                }
            }
        }
        PhysicalContactGenerator {
            dat,
            contact_bbs,
            contacts,
            params,
        }
    }

    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn close_edges_on_cell(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        cell_polys: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let v = cell_polys[ci][vi];
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| {
                let (t, d) = self.dat.get(ci, vi, oci, ovi);
                if d < f32::INFINITY {
                    let p0 = cell_polys[oci][ovi];
                    let p1 =
                        cell_polys[oci][circ_ix_plus(ovi, NVERTS)];
                    let p = t * (p1 - p0) + p0;
                    let delta = p - v;

                    let edge_rgtp = (cell_rgtps[oci][ovi]
                        + cell_rgtps[oci][circ_ix_plus(ovi, NVERTS)])
                        / 2.0;
                    Some(CloseEdge {
                        cell_ix: oci,
                        vert_ix: ovi,
                        crl: CrlEffect::calc(v_rgtp, edge_rgtp),
                        delta,
                        t,
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<CloseEdge>>()
    }

    /// Get edges which contain points close to vertex `vi` on cell `ci`.
    pub fn close_edges(
        &self,
        ci: usize,
        vi: usize,
        cell_polys: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.dat.num_cells {
            r.append(&mut self.close_edges_on_cell(
                ci, vi, oci, cell_polys, cell_rgtps,
            ))
        }
        r
    }

    pub fn update(
        &mut self,
        ci: usize,
        bb: &BBox,
        vs: &[V2D; NVERTS],
        cell_polys: &[[V2D; NVERTS]],
    ) {
        let bb = bb.expand_by(self.params.range);
        self.contact_bbs[ci] = bb;
        // for (oxi, obb) in
        //     self.contact_bbs[(ci + 1)..].iter().enumerate()
        // {
        //     self.contacts.set(
        //         ci,
        //         ci + 1 + oxi,
        //         obb.intersects(&bb),
        //     );
        // }
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contacts.set(ci, oci, obb.intersects(&bb));
            }
        }
        for (oci, ovs) in cell_polys.iter().enumerate() {
            if ci != oci && self.contacts.get(ci, oci) {
                for (vi, v) in vs.iter().enumerate() {
                    for (ovi, ov) in ovs.iter().enumerate() {
                        let ov2 = &ovs[circ_ix_plus(ovi, NVERTS)];
                        let (t, d) = seg_to_point_dist(v, ov, ov2);
                        match d.partial_cmp(&self.params.range) {
                            Some(ord) => match ord {
                                Ordering::Greater => {
                                    continue;
                                }
                                _ => {
                                    self.dat.set(
                                        ci,
                                        vi,
                                        oci,
                                        ovi,
                                        (t, d),
                                    );
                                    break;
                                }
                            },
                            _ => panic!(
                                "cannot compare {} and  {}",
                                d, self.params.range
                            ),
                        }
                    }
                }
            }
        }
    }

    pub fn generate(
        &self,
        cell_polys: &[[V2D; NVERTS]],
        all_rgtps: &[[f32; NVERTS]],
    ) -> PhysContactFactors {
        let num_cells = self.contacts.num_cells;
        let mut adh = vec![[V2D::default(); NVERTS]; num_cells];
        let mut cal = vec![[0.0f32; NVERTS]; num_cells];
        let mut cil = vec![[0.0f32; NVERTS]; num_cells];
        for ci in 0..num_cells {
            let x_cals = &mut cal[ci];
            let x_cils = &mut cil[ci];
            for vi in 0..NVERTS {
                for oci in 0..num_cells {
                    if oci != ci {
                        for CloseEdge {
                            cell_ix: oci,
                            vert_ix: ovi,
                            crl,
                            delta,
                            t,
                        } in self
                            .close_edges(
                                ci, vi, cell_polys, all_rgtps,
                            )
                            .into_iter()
                        {
                            match (self.params.cal_mag, crl) {
                                (Some(cal_mag), CrlEffect::Cal) => {
                                    x_cals[vi] = cal_mag;
                                }
                                (Some(_), CrlEffect::Cil)
                                | (None, _) => {
                                    x_cils[vi] = self.params.cil_mag;
                                }
                            }

                            match self.params.adh_mag {
                                Some(adh_mag) => {
                                    let d = delta.mag();
                                    if d.abs() > 1e-3
                                        && d < self.params.range
                                    {
                                        let adh_force = adh_mag
                                            * d
                                            * delta.unitize();
                                        adh[oci][ovi] = adh[oci][ovi]
                                            + -1.0
                                                * (1.0 - t)
                                                * adh_force;
                                        adh[oci][circ_ix_plus(
                                            ovi, NVERTS,
                                        )] =
                                            adh[oci][circ_ix_plus(
                                                ovi, NVERTS,
                                            )] + -1.0 * t * adh_force;
                                        adh[ci][vi] =
                                            adh[ci][vi] + adh_force;
                                    }
                                }
                                None => {}
                            };
                        }
                    }
                }
            }
        }
        PhysContactFactors { adh, cil, cal }
    }
}

pub enum CrlEffect {
    Cil,
    Cal,
}

impl CrlEffect {
    pub fn calc(a: RgtpState, b: RgtpState) -> CrlEffect {
        match (a > 0.0, b > 0.0) {
            (true, true) => CrlEffect::Cal,
            (_, _) => CrlEffect::Cil,
        }
    }
}

/// Information relevant to calculating CIL and adhesion due to an edge in proximity to the focus
/// vertex.
pub struct CloseEdge {
    /// Cell containing edge close to focus vertex.
    pub cell_ix: usize,
    /// Close edge runs from `vert_ix` to `vert_ix + 1`.
    pub vert_ix: usize,
    /// Contact regulation of motion.
    pub crl: CrlEffect,
    /// Let the position of the point on the close edge closest to the focus vertex be denoted `p`,
    /// and the position of the focus vertex be denoted `v`. `delta` is such that `delta + v = p`.
    pub delta: V2D,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub t: f32,
}
