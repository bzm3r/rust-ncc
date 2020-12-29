use crate::interactions::dat4d::CvCvDat;
use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::{
    generate_contacts, ContactData, RgtpState,
};
use crate::math::close_to_zero;
use crate::math::geometry::BBox;
use crate::math::v2d::V2D;
use crate::parameters::PhysicalContactParams;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use std::cmp::Ordering;

#[derive(Clone, Copy)]
pub struct Dist(f32);
#[derive(Clone, Copy)]
pub struct LineSegParam(f32);

#[derive(Clone, Copy)]
pub enum ClosePoint {
    Vertex(V2D),
    OnEdge(LineSegParam, V2D),
    None,
}

impl ClosePoint {
    /// Returns the point closest to `p` (`ClosePoint`) on the line
    /// segment `k = (b - a)*t + a, 0 <= t < 1`.
    pub fn calc(range: f32, p: V2D, a: V2D, b: V2D) -> ClosePoint {
        // Is `p` close to `a`? Then it interacts directly with `a`.
        let pa = a - p;
        if pa.mag() < range {
            ClosePoint::Vertex(pa)
        } else if (b - p).mag() < range {
            ClosePoint::None
        } else {
            let seg_gen = b - a;
            let t = seg_gen.dot(&p) / seg_gen.mag_squared();
            // Is `t` in the interval `[0, 1)`? If yes, then the close
            // point lies on the edge.
            if {
                match t.partial_cmp(&0.0) {
                    Some(Ordering::Less) => false,
                    Some(Ordering::Greater)
                    | Some(Ordering::Equal) => true,
                    _ => panic!(
                        "partial_cmp could not determine ordering of t ({}) relative to 0.0", t
                    ),
                }
            } && (t < 1.0)
            {
                let pc = t * (seg_gen + a) - p;
                if pc.mag() < range {
                    ClosePoint::OnEdge(LineSegParam(t), pc)
                } else {
                    ClosePoint::None
                }
            } else {
                ClosePoint::None
            }
        }
    }
}

/// Generates CIL/CAL/adhesion related interaction information. These
/// are the interactions that require cells to engage in
/// physical contact.
#[derive(Clone)]
pub struct PhysicalContactGenerator {
    dat: CvCvDat<ClosePoint>,
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
        let mut dat = CvCvDat::empty(num_cells, ClosePoint::None);
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
                            let ou = &ovs[circ_ix_plus(ovi, NVERTS)];
                            dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                ClosePoint::calc(
                                    params.range,
                                    *v,
                                    *ov,
                                    *ou,
                                ),
                            )
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
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| match self.dat.get(ci, vi, oci, ovi) {
                ClosePoint::None => None,
                ClosePoint::OnEdge(LineSegParam(t), delta) => {
                    //TODO: confirm that we don't want:
                    // let edge_rgtp = (1.0 - t) * cell_rgtps[oci][ovi]
                    //     + t * cell_rgtps[oci]
                    //         [circ_ix_plus(ovi, NVERTS)];
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
                }
                ClosePoint::Vertex(delta) => Some(CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl: CrlEffect::calc(
                        v_rgtp,
                        cell_rgtps[oci][ovi],
                    ),
                    delta,
                    t: 0.0,
                }),
            })
            .collect::<Vec<CloseEdge>>()
    }

    /// Get edges which contain points close to vertex `vi` on cell `ci`.
    pub fn close_edges(
        &self,
        ci: usize,
        vi: usize,
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.dat.num_cells {
            r.append(
                &mut self
                    .close_edges_on_cell(ci, vi, oci, cell_rgtps),
            )
        }
        r
    }

    /// Get vertices on cell `oci` that are close to cell `ci`.
    pub fn get_close_verts(
        &self,
        ci: usize,
        oci: usize,
    ) -> Vec<usize> {
        (0..NVERTS).collect()
        // let mut r = vec![];
        // for vi in 0..NVERTS {
        //     for (ovi, cp) in self
        //         .dat
        //         .get_per_other_vertex(ci, vi, oci)
        //         .iter()
        //         .enumerate()
        //     {
        //         match cp {
        //             ClosePoint::Vertex(_)
        //             | ClosePoint::OnEdge(_, _) => r.push(ovi),
        //             ClosePoint::None => {}
        //         }
        //     }
        // }
        // r.dedup();
        // r.sort_unstable();
        // r.dedup();
        // r
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
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contacts.set(ci, oci, obb.intersects(&bb));
            }
        }
        for (oci, ovs) in cell_polys.iter().enumerate() {
            if ci != oci && self.contacts.get(ci, oci) {
                for (vi, v) in vs.iter().enumerate() {
                    for (ovi, ov) in ovs.iter().enumerate() {
                        let ow = &ovs[circ_ix_plus(ovi, NVERTS)];
                        self.dat.set(
                            ci,
                            vi,
                            oci,
                            ovi,
                            ClosePoint::calc(
                                self.params.range,
                                *v,
                                *ov,
                                *ow,
                            ),
                        );
                    }
                }
            }
        }
    }

    pub fn generate(
        &self,
        cell_rgtps: &[[f32; NVERTS]],
    ) -> PhysContactFactors {
        let num_cells = self.contacts.num_cells;
        let mut adh = vec![[V2D::default(); NVERTS]; num_cells];
        let mut cal = vec![[0.0f32; NVERTS]; num_cells];
        let mut cil = vec![[0.0f32; NVERTS]; num_cells];
        for ci in 0..num_cells {
            let x_cals = &mut cal[ci];
            let x_cils = &mut cil[ci];
            for vi in 0..NVERTS {
                for CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl,
                    delta,
                    t,
                    ..
                } in
                    self.close_edges(ci, vi, cell_rgtps).into_iter()
                {
                    match (self.params.cal_mag, crl) {
                        (Some(cal_mag), CrlEffect::Cal) => {
                            x_cals[vi] = cal_mag;
                        }
                        (Some(_), CrlEffect::Cil) | (None, _) => {
                            x_cils[vi] = self.params.cil_mag;
                        }
                    }

                    if let Some(adh_mag) = self.params.adh_mag {
                        let adh_force = adh_mag
                            * ((1.0 / self.params.range) * delta);
                        //* ((1.0 / self.params.range) * delta);
                        // We are close to the vertex.
                        if close_to_zero(t) {
                            adh[oci][ovi] = adh[oci][ovi] - adh_force;
                        } else {
                            adh[oci][ovi] = adh[oci][ovi]
                                + -1.0 * (1.0 - t) * adh_force;
                            let owi = circ_ix_plus(ovi, NVERTS);
                            adh[oci][owi] =
                                adh[oci][owi] + -1.0 * t * adh_force;
                            adh[ci][vi] = adh[ci][vi] + adh_force;
                        }
                    };
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
