use crate::interactions::dat_4d::CvCvDat;
use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::{generate_contacts, RgtpActivityDiff};
use crate::math::geometry::{BBox, Poly};
use crate::math::v2d::V2D;
use crate::math::{
    capped_linear_fn, close_to_zero, in_unit_interval,
};
use crate::parameters::{CloseBounds, PhysicalContactParams};
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy)]
pub struct Dist(f32);
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct LineSegParam(f32);

#[derive(Clone, Copy, Serialize, Deserialize)]
pub enum ClosePoint {
    Vertex {
        vector_to: V2D,
        smooth_factor: f32,
    },
    OnEdge {
        edge_point_param: f32,
        vector_to: V2D,
        smooth_factor: f32,
    },
    None,
}

impl ClosePoint {
    /// Returns the point closest to `p` (`ClosePoint`) on the line
    /// segment `k = (b - a)*t + a, 0 <= t < 1`.
    pub fn calc(
        range: CloseBounds,
        p: V2D,
        a: V2D,
        b: V2D,
    ) -> ClosePoint {
        // Is `p` close to `a`? Then it interacts directly with `a`.
        let pa = a - p;
        if pa.mag() < range.zero_at {
            let smooth_factor = capped_linear_fn(
                pa.mag(),
                range.zero_at,
                range.one_at,
            );
            ClosePoint::Vertex {
                vector_to: pa,
                smooth_factor,
            }
        } else if (b - p).mag() < range.zero_at {
            ClosePoint::None
        } else {
            let ab = b - a;
            let t = ab.dot(&pa) / ab.mag();
            // Is `t` in the interval `[0, 1)`? If yes, then the close
            // point lies on the edge.
            if in_unit_interval(t) {
                let pc = (t * ab) + pa;
                if pc.mag() < range.zero_at {
                    ClosePoint::OnEdge {
                        edge_point_param: t,
                        vector_to: pc,
                        smooth_factor: capped_linear_fn(
                            pc.mag(),
                            range.zero_at,
                            range.one_at,
                        ),
                    }
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
#[derive(Clone, Deserialize, Serialize)]
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
        cell_polys: &[Poly],
        params: PhysicalContactParams,
    ) -> PhysicalContactGenerator {
        let num_cells = cell_polys.len();
        let mut dat = CvCvDat::empty(num_cells, ClosePoint::None);
        let contact_bbs = cell_polys
            .iter()
            .map(|cp| cp.bbox.expand_by(params.range.zero_at))
            .collect::<Vec<BBox>>();
        let contacts = generate_contacts(&contact_bbs);
        let c0_contacts = contacts.get(0, 1);
        let c1_contacts = contacts.get(1, 0);
        for (ai, poly) in cell_polys.iter().enumerate() {
            for (bi, other) in cell_polys.iter().enumerate() {
                if ai != bi && contacts.get(ai, bi) {
                    for (avi, p) in poly.verts.iter().enumerate() {
                        for (bvi, a) in other.verts.iter().enumerate()
                        {
                            let b = &other.verts
                                [circ_ix_plus(bvi, NVERTS)];
                            let cp = ClosePoint::calc(
                                params.range,
                                *p,
                                *a,
                                *b,
                            );
                            dat.set(
                                ai,
                                avi,
                                bi,
                                bvi,
                                ClosePoint::calc(
                                    params.range,
                                    *p,
                                    *a,
                                    *b,
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
        cell_rgtps: &[[RgtpActivityDiff; NVERTS]],
    ) -> Vec<CloseEdge> {
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| match self.dat.get(ci, vi, oci, ovi) {
                ClosePoint::None => None,
                ClosePoint::OnEdge {
                    vector_to,
                    smooth_factor,
                    edge_point_param,
                } => {
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
                        vector_to,
                        edge_point_param,
                        smooth_factor,
                    })
                }
                ClosePoint::Vertex {
                    vector_to,
                    smooth_factor,
                } => Some(CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl: CrlEffect::calc(
                        v_rgtp,
                        cell_rgtps[oci][ovi],
                    ),
                    vector_to,
                    edge_point_param: 0.0,
                    smooth_factor,
                }),
            })
            .collect::<Vec<CloseEdge>>()
    }

    /// Get edges which contain points close to vertex `vi` on cell `ci`.
    pub fn close_edges(
        &self,
        ci: usize,
        vi: usize,
        cell_rgtps: &[[RgtpActivityDiff; NVERTS]],
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

    // /// Get vertices on cell `oci` that are close to cell `ci`.
    // pub fn get_close_verts(
    //     &self,
    //     aci: usize,
    //     bci: usize,
    // ) -> Vec<usize> {
    //     let mut r = vec![];
    //     for avi in 0..NVERTS {
    //         for (bvi, close_point) in self
    //             .dat
    //             .get_per_other_vertex(aci, avi, bci)
    //             .iter()
    //             .enumerate()
    //         {
    //             match close_point {
    //                 ClosePoint::Vertex(_)
    //                 | ClosePoint::OnEdge(_, _) => r.push(bvi),
    //                 ClosePoint::None => {}
    //             }
    //         }
    //     }
    //     r.sort_unstable();
    //     r.dedup();
    //     r
    // }

    pub fn update(&mut self, ci: usize, cell_polys: &[Poly]) {
        let poly = cell_polys[ci];
        let bb =
            cell_polys[ci].bbox.expand_by(self.params.range.zero_at);
        self.contact_bbs[ci] = bb;
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contacts.set(ci, oci, obb.intersects(&bb));
            }
        }
        for (oci, other) in cell_polys.iter().enumerate() {
            if ci != oci && self.contacts.get(ci, oci) {
                for (pi, p) in poly.verts.iter().enumerate() {
                    for (ai, a) in other.verts.iter().enumerate() {
                        if ci == 0 && oci == 1 && pi == 3 && ai == 13
                        {
                            // println!(
                            //     "c0_3, c1_13: {}",
                            //     (p - a).mag()
                            // );
                        }
                        if ci == 0 && oci == 1 && pi == 5 && ai == 11
                        {
                            // println!(
                            //     "c0_5, c1_11: {}",
                            //     (p - a).mag()
                            // );
                        }
                        if ci == 0 && oci == 1 && pi == 4 && ai == 12
                        {
                            println!(
                                "c0_4, c1_12: {}",
                                (p - a).mag()
                            );
                        }
                        if ci == 1 && oci == 0 && pi == 12 && ai == 4
                        {
                            println!(
                                "c1_12, c0_4: {}",
                                (p - a).mag()
                            );
                        }
                        let bi = circ_ix_plus(ai, NVERTS);
                        let b = &other.verts[bi];
                        self.dat.set(
                            ci,
                            pi,
                            oci,
                            ai,
                            ClosePoint::calc(
                                self.params.range,
                                *p,
                                *a,
                                *b,
                            ),
                        );
                    }
                }
            }
        }
    }

    pub fn generate(
        &self,
        cell_rgtps: &[[RgtpActivityDiff; NVERTS]],
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
                    vector_to,
                    edge_point_param,
                    smooth_factor,
                } in
                    self.close_edges(ci, vi, cell_rgtps).into_iter()
                {
                    match (self.params.cal_mag, crl) {
                        (Some(cal_mag), CrlEffect::Cal) => {
                            x_cals[vi] = smooth_factor * cal_mag;
                        }
                        (Some(_), CrlEffect::Cil) | (None, _) => {
                            x_cils[vi] =
                                smooth_factor * self.params.cil_mag;
                        }
                    }

                    if let Some(adh_mag) = self.params.adh_mag {
                        let d = vector_to.mag();
                        let s = capped_linear_fn(
                            d,
                            self.params.range.one_at,
                            0.0,
                        );
                        let adh_force = adh_mag * s * vector_to;
                        //* ((1.0 / self.params.range) * delta);
                        // We are close to the vertex.
                        if close_to_zero(edge_point_param) {
                            adh[oci][ovi] = adh[oci][ovi] - adh_force;
                        } else {
                            adh[oci][ovi] = adh[oci][ovi]
                                + -1.0
                                    * (1.0 - edge_point_param)
                                    * adh_force;
                            let owi = circ_ix_plus(ovi, NVERTS);
                            adh[oci][owi] = adh[oci][owi]
                                + -1.0 * edge_point_param * adh_force;
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
    pub fn calc(
        a: RgtpActivityDiff,
        b: RgtpActivityDiff,
    ) -> CrlEffect {
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
    /// Let the position of the focus vertex be denoted `p`, and
    /// the point on the close edge closest to the focus vertex
    /// be denoted `c`. `delta` is such that `delta + p = c`.
    pub vector_to: V2D,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub edge_point_param: f32,
    pub smooth_factor: f32,
}
