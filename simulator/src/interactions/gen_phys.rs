use crate::interactions::dat_4d::CvCvDat;
use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::{gen_contact_matrix, RelativeRgtpActivity};
use crate::math::geometry::{BBox, Poly};
use crate::math::v2d::V2d;
use crate::math::{
    capped_linear_fn, close_to_zero, in_unit_interval, InUnitInterval,
};
use crate::parameters::PhysicalContactParams;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Copy)]
pub struct Dist(f64);
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct LineSegParam(f64);

#[derive(Clone, Copy, Serialize, Deserialize, Debug, PartialEq)]
pub enum ClosePoint {
    Vertex {
        vector_to: V2d,
        smooth_factor: f64,
        dist_sq: f64,
    },
    OnEdge {
        edge_point_param: f64,
        vector_to: V2d,
        smooth_factor: f64,
        dist_sq: f64,
    },
    None {
        dist_sq: f64,
    },
}

impl Default for ClosePoint {
    fn default() -> Self {
        ClosePoint::None {
            dist_sq: f64::INFINITY,
        }
    }
}

impl ClosePoint {
    /// Returns the point closest to `p` (`ClosePoint`) on the line
    /// segment `k = (b - a)*t + a, 0 <= t < 1`.
    pub fn calc(
        one_at: f64,
        zero_at: f64,
        zero_at_sq: f64,
        test_point: V2d,
        seg_start: V2d,
        seg_end: V2d,
    ) -> ClosePoint {
        // Is `p` close to `a`? Then it interacts directly with `a`.
        let s_to_tp = test_point - seg_start;
        let s_to_tp_dsq = s_to_tp.mag_squared();
        if s_to_tp_dsq < 1e-6 {
            let smooth_factor =
                capped_linear_fn(s_to_tp_dsq.sqrt(), zero_at, one_at);
            ClosePoint::Vertex {
                vector_to: -1.0 * s_to_tp,
                smooth_factor,
                dist_sq: s_to_tp_dsq,
            }
        } else {
            let e_to_tp_dsq = (seg_end - test_point).mag_squared();
            if e_to_tp_dsq < 1e-6 {
                ClosePoint::None {
                    dist_sq: e_to_tp_dsq,
                }
            } else {
                let seg_vec = seg_end - seg_start;
                let t = seg_vec.dot(&s_to_tp) / seg_vec.mag_squared();
                // Is `t` in the interval `[0, 1)`? If yes, then the close
                // point lies on the edge.
                match in_unit_interval(t, 1e-3) {
                    InUnitInterval::Zero | InUnitInterval::In => {
                        let c = t * (seg_vec) + seg_start;
                        let tp_to_c = c - test_point;
                        let tp_to_c_mag_sq = tp_to_c.mag_squared();
                        if tp_to_c_mag_sq < zero_at_sq {
                            ClosePoint::OnEdge {
                                edge_point_param: t,
                                vector_to: tp_to_c,
                                smooth_factor: capped_linear_fn(
                                    tp_to_c_mag_sq.sqrt(),
                                    zero_at,
                                    one_at,
                                ),
                                dist_sq: tp_to_c_mag_sq,
                            }
                        } else {
                            ClosePoint::None {
                                dist_sq: tp_to_c_mag_sq,
                            }
                        }
                    }
                    _ => ClosePoint::None {
                        dist_sq: f64::INFINITY,
                    },
                }
            }
        }
    }
}

impl fmt::Display for ClosePoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClosePoint::Vertex {
                vector_to,
                smooth_factor,
                dist_sq,
            } => {
                write!(
                    f,
                    "Vertex(vector_to: {}, smooth_factor: {}, dist: {})",
                    vector_to.mag(),
                    smooth_factor,
                    dist_sq.sqrt(),
                )
            }
            ClosePoint::OnEdge {
                edge_point_param,
                vector_to,
                smooth_factor,
                dist_sq,
            } => {
                write!(f, "OnEdge(edge_point_param: {}, vector_to: {}, smooth_factor: {}, \
                dist: {})", edge_point_param, vector_to.mag(), smooth_factor, dist_sq.sqrt())
            }
            ClosePoint::None { dist_sq } => {
                write!(f, "None(dist: {})", dist_sq.sqrt())
            }
        }
    }
}

/// Generates CIL/CAL/adhesion related interaction information. These
/// are the interactions that require cells to engage in
/// physical contact.
#[derive(
    Clone, Deserialize, Serialize, PartialEq, Default, Debug,
)]
pub struct PhysicalContactGenerator {
    dat: CvCvDat<ClosePoint>,
    min_dist_matrix: Vec<[f64; NVERTS]>,
    pub contact_bbs: Vec<BBox>,
    pub contact_matrix: SymCcDat<bool>,
    pub params: PhysicalContactParams,
}

pub struct PhysContactFactors {
    pub adh: Vec<[V2d; NVERTS]>,
    pub cil: Vec<[f64; NVERTS]>,
    pub cal: Vec<[f64; NVERTS]>,
}

impl PhysicalContactGenerator {
    /// Calculate distances between vertices of cells in contact.
    pub fn new(
        cell_polys: &[Poly],
        params: PhysicalContactParams,
    ) -> PhysicalContactGenerator {
        let num_cells = cell_polys.len();
        let mut dat =
            CvCvDat::empty(num_cells, ClosePoint::default());
        let contact_bbs = cell_polys
            .iter()
            .map(|cp| cp.bbox.expand_by(params.zero_at_sq))
            .collect::<Vec<BBox>>();
        let contact_matrix = gen_contact_matrix(&contact_bbs);
        for (ci, poly) in cell_polys.iter().enumerate() {
            for (oci, other) in cell_polys.iter().enumerate() {
                if ci != oci {
                    for (vi, v) in poly.verts.iter().enumerate() {
                        for (ovi, ov) in
                            other.verts.iter().enumerate()
                        {
                            let ow = &other.verts
                                [circ_ix_plus(ovi, NVERTS)];
                            if contact_matrix.get(ci, oci) {
                                dat.set(
                                    ci,
                                    vi,
                                    oci,
                                    ovi,
                                    ClosePoint::calc(
                                        params.one_at,
                                        params.zero_at,
                                        params.zero_at_sq,
                                        *v,
                                        *ov,
                                        *ow,
                                    ),
                                );
                            } else {
                                dat.set(
                                    ci,
                                    vi,
                                    oci,
                                    ovi,
                                    ClosePoint::default(),
                                )
                            }
                        }
                    }
                }
            }
        }
        let mut min_dist_matrix =
            vec![[f64::INFINITY; NVERTS]; num_cells];
        for ci in 0..num_cells {
            min_dist_matrix[ci] =
                PhysicalContactGenerator::eval_min_dist(
                    ci,
                    &contact_matrix,
                    &dat,
                );
        }
        PhysicalContactGenerator {
            dat,
            min_dist_matrix,
            contact_bbs,
            contact_matrix,
            params,
        }
    }

    pub fn update(&mut self, ci: usize, cell_polys: &[Poly]) {
        let poly = cell_polys[ci];
        let bb =
            cell_polys[ci].bbox.expand_by(self.params.zero_at_sq);
        self.contact_bbs[ci] = bb;
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != ci {
                self.contact_matrix.set(ci, oci, obb.intersects(&bb));
            }
        }
        for (oci, other) in cell_polys.iter().enumerate() {
            if ci != oci {
                for (vi, v) in poly.verts.iter().enumerate() {
                    let w = &poly.verts[circ_ix_plus(vi, NVERTS)];
                    for (ovi, ov) in other.verts.iter().enumerate() {
                        let ow =
                            &other.verts[circ_ix_plus(ovi, NVERTS)];
                        if self.contact_matrix.get(ci, oci) {
                            self.dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                ClosePoint::calc(
                                    self.params.one_at,
                                    self.params.zero_at,
                                    self.params.zero_at_sq,
                                    *v,
                                    *ov,
                                    *ow,
                                ),
                            );
                            self.dat.set(
                                oci,
                                ovi,
                                ci,
                                vi,
                                ClosePoint::calc(
                                    self.params.one_at,
                                    self.params.zero_at,
                                    self.params.zero_at_sq,
                                    *ov,
                                    *v,
                                    *w,
                                ),
                            );
                        } else {
                            self.dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                ClosePoint::default(),
                            );
                            self.dat.set(
                                oci,
                                ovi,
                                ci,
                                vi,
                                ClosePoint::default(),
                            );
                        }
                    }
                }
            }
        }
        self.min_dist_matrix[ci] =
            PhysicalContactGenerator::eval_min_dist(
                ci,
                &self.contact_matrix,
                &self.dat,
            );
    }

    pub fn eval_min_dist(
        ci: usize,
        contact_matrix: &SymCcDat<bool>,
        dat: &CvCvDat<ClosePoint>,
    ) -> [f64; NVERTS] {
        let mut r = [f64::INFINITY; NVERTS];
        for vi in 0..NVERTS {
            for oci in 0..contact_matrix.num_cells {
                if oci != ci && contact_matrix.get(ci, oci) {
                    for ovi in 0..NVERTS {
                        match dat.get(ci, vi, oci, ovi) {
                            ClosePoint::OnEdge {
                                dist_sq, ..
                            }
                            | ClosePoint::Vertex {
                                dist_sq, ..
                            }
                            | ClosePoint::None { dist_sq } => {
                                r[vi] = r[vi].min(dist_sq);
                            }
                        }
                    }
                }
            }
        }
        r
    }

    pub fn min_dist_to(&self, ci: usize, vi: usize) -> f64 {
        self.min_dist_matrix[ci][vi]
    }

    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn get_close_edges_on_cell(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        rel_rgtps_per_cell: &[[RelativeRgtpActivity; NVERTS]],
    ) -> Vec<CloseEdge> {
        let v_rgtp = rel_rgtps_per_cell[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| match self.dat.get(ci, vi, oci, ovi) {
                ClosePoint::None { .. } => None,
                ClosePoint::OnEdge {
                    vector_to,
                    smooth_factor,
                    edge_point_param,
                    ..
                } => {
                    let edge_rgtp = RelativeRgtpActivity::mix_rel_rgtp_act_across_edge(
                        rel_rgtps_per_cell[oci][ovi],
                        rel_rgtps_per_cell[oci]
                            [circ_ix_plus(ovi, NVERTS)], edge_point_param,
                    );
                    Some(CloseEdge {
                        cell_ix: oci,
                        vert_ix: ovi,
                        crl: CrlEffect::calc_crl_on_focus(
                            v_rgtp, edge_rgtp,
                        ),
                        vector_to,
                        edge_point_param,
                        smooth_factor,
                    })
                }
                ClosePoint::Vertex {
                    vector_to,
                    smooth_factor,
                    ..
                } => Some(CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl: CrlEffect::calc_crl_on_focus(
                        v_rgtp,
                        rel_rgtps_per_cell[oci][ovi],
                    ),
                    vector_to,
                    edge_point_param: 0.0,
                    smooth_factor,
                }),
            })
            .collect::<Vec<CloseEdge>>()
    }

    /// Get edges which contain points close to vertex `vi` on cell `ci`.
    pub fn get_close_edges_to(
        &self,
        ci: usize,
        vi: usize,
        cell_rgtps: &[[RelativeRgtpActivity; NVERTS]],
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.dat.num_cells {
            r.append(
                &mut self
                    .get_close_edges_on_cell(ci, vi, oci, cell_rgtps),
            )
        }
        r
    }

    pub fn generate(
        &self,
        rel_rgtps_per_cell: &[[RelativeRgtpActivity; NVERTS]],
    ) -> PhysContactFactors {
        let num_cells = self.contact_matrix.num_cells;
        let mut adh_per_cell =
            vec![[V2d::default(); NVERTS]; num_cells];
        let mut cal_per_cell = vec![[0.0f64; NVERTS]; num_cells];
        let mut cil_per_cell = vec![[0.0f64; NVERTS]; num_cells];
        for ci in 0..num_cells {
            let x_cals = &mut cal_per_cell[ci];
            let x_cils = &mut cil_per_cell[ci];
            for vi in 0..NVERTS {
                for CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl,
                    vector_to,
                    edge_point_param,
                    smooth_factor,
                } in self
                    .get_close_edges_to(ci, vi, rel_rgtps_per_cell)
                    .into_iter()
                {
                    match (self.params.cal_mag, crl) {
                        (Some(cal_mag), CrlEffect::Cal) => {
                            x_cals[vi] = x_cals[vi]
                                .max(smooth_factor * cal_mag);
                        }
                        (Some(_), CrlEffect::Cil) | (None, _) => {
                            x_cils[vi] = x_cils[vi].max(
                                smooth_factor * self.params.cil_mag,
                            );
                        }
                    }

                    if let Some(adh_mag) = self.params.adh_mag {
                        let vc_mag = vector_to.mag();
                        let adh_strain = if vc_mag
                            > self.params.adh_max
                        {
                            if vc_mag > self.params.zero_at {
                                0.0
                            } else {
                                1.0 - ((vc_mag - self.params.adh_max)
                                    / self.params.adh_delta_break)
                            }
                        } else {
                            (vc_mag / self.params.adh_rest) - 1.0
                        };
                        if ci == 0 && vi == 0 && oci == 1 && ovi == 8
                        {
                            println!("zero_at: {}, adh_max: {}, adh_delta_break: {}, adh_rest: {}, vc_mag: {}, adh_strain: {}", self.params
                                .zero_at, self.params.adh_max, self.params.adh_delta_break, self
                                         .params.adh_rest, vc_mag, adh_strain);
                        }
                        let adh_force = adh_mag
                            * adh_strain
                            * vector_to.unitize();
                        //* ((1.0 / self.params.range) * delta);
                        // We are close to the vertex.
                        if close_to_zero(edge_point_param, 1e-3) {
                            adh_per_cell[oci][ovi] =
                                adh_per_cell[oci][ovi] - adh_force;
                            adh_per_cell[ci][vi] =
                                adh_per_cell[ci][vi] + adh_force;
                        } else {
                            adh_per_cell[oci][ovi] = adh_per_cell
                                [oci][ovi]
                                - (1.0 - edge_point_param)
                                    * adh_force;
                            let owi = circ_ix_plus(ovi, NVERTS);
                            adh_per_cell[oci][owi] = adh_per_cell
                                [oci][owi]
                                - edge_point_param * adh_force;
                            adh_per_cell[ci][vi] =
                                adh_per_cell[ci][vi] + adh_force;
                        }
                    };
                }
            }
        }
        PhysContactFactors {
            adh: adh_per_cell,
            cil: cil_per_cell,
            cal: cal_per_cell,
        }
    }
}

pub enum CrlEffect {
    Cil,
    Cal,
}

impl CrlEffect {
    //TODO: should CIL/CAL be modelled with a "relative strength"?
    pub fn calc_crl_on_focus(
        focus_vertex: RelativeRgtpActivity,
        other: RelativeRgtpActivity,
    ) -> CrlEffect {
        use RelativeRgtpActivity::{RacDominant, RhoDominant};
        match (focus_vertex, other) {
            (RacDominant(_), RhoDominant(_)) => CrlEffect::Cal,
            (RhoDominant(_), RacDominant(_))
            | (RhoDominant(_), RhoDominant(_))
            | (RacDominant(_), RacDominant(_)) => CrlEffect::Cil,
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
    pub vector_to: V2d,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub edge_point_param: f64,
    pub smooth_factor: f64,
}
