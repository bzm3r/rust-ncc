use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::dat_sym4d::SymCcVvDat;
use crate::interactions::gen_contact_matrix;
use crate::interactions::gen_phys::PhysicalContactGenerator;
use crate::math::geometry::{BBox, LineSeg2D, Poly};
use crate::parameters::CoaParams;
use crate::utils::circ_ix_minus;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Deserialize, Serialize)]
pub struct VertexPairInfo {
    dist: f64,
    num_intersects: f64,
}

impl VertexPairInfo {
    pub fn infinity() -> VertexPairInfo {
        VertexPairInfo {
            dist: f64::INFINITY,
            num_intersects: f64::INFINITY,
        }
    }
}

#[derive(Clone, Deserialize, Serialize)]
pub struct CoaGenerator {
    dat: SymCcVvDat<VertexPairInfo>,
    contact_bbs: Vec<BBox>,
    contact_matrix: SymCcDat<bool>,
    cell_polys: Vec<Poly>,
    coa_matrix: Vec<[f64; NVERTS]>,
    params: CoaParams,
}

/// Suppose we have a line segment `lseg`, which goes from vertex `vi_a` of cell `A`
/// to vertex `vi_b` of cell `B`. `lseg` models a COA interaction between cell `A`
/// and `B` due to interaction between vertex `vi_a` and `vi_b`.
///
/// Let `C` be a cell that is not `A` or `B` (an "other" polygon). This function
/// checks to see if `lseg` intersects `C`. Note that `lseg` can only intersect `C` in
/// the following ways: 1. `lseg` passes through two of the vertices of `C` (weak
/// intersection). 2. `lseg` intersects an edge of `C` at a point that is not also
/// a vertex of `C` (strong intersection).
///
/// We are restricted to these cases because `C` is by assumption a cell which
/// does not contain the endpoints of `lseg`. Note that case `1` is unlikely, but
/// possible especially in the initial if cells have been initialized in a
/// regular lattice.
pub fn check_other_poly_intersect(
    lseg: &LineSeg2D,
    poly: &Poly,
) -> bool {
    if lseg.intersects_bbox(&poly.bbox) {
        for edge in poly.edges.iter() {
            if lseg.check_intersection(edge) {
                return true;
            }
        }
    }
    false
}

/// Suppose we have a line segment `lseg`, which goes from vertex `vi_a` of cell `A`
/// to vertex `vi_b` of cell `B`. `lseg` models a COA interaction between cell `A`
/// and `B` due to interaction between vertex `vi_a` and `vi_b`.
///
/// It could be that that `lseg` one of `A` or `B`, the "root" polygons of `lseg`.
/// This function checks if this has occurred, but ignores intersections
/// involving the source/destination vertices of `lseg`.
pub fn check_root_poly_intersect(
    lseg: &LineSeg2D,
    poly_a: &Poly,
    poly_b: &Poly,
    vi_a: usize,
    vi_b: usize,
) -> bool {
    let ui_a = circ_ix_minus(vi_a, NVERTS);
    for (ei, edge) in poly_a.edges.iter().enumerate() {
        if ei != vi_a && ei != ui_a && edge.check_intersection(lseg) {
            return true;
        }
    }

    let ui_b = circ_ix_minus(vi_b, NVERTS);
    for (ei, edge) in poly_b.edges.iter().enumerate() {
        if ei != vi_b && ei != ui_b && edge.check_intersection(lseg) {
            return true;
        }
    }

    false
}

/// Calculate clearance and distance.
pub fn calc_pair_info(
    ci: usize,
    vi: usize,
    oci: usize,
    ovi: usize,
    lseg: LineSeg2D,
    cell_polys: &[Poly],
) -> VertexPairInfo {
    if check_root_poly_intersect(
        &lseg,
        &cell_polys[ci],
        &cell_polys[oci],
        vi,
        ovi,
    ) {
        return VertexPairInfo {
            dist: lseg.len,
            num_intersects: f64::INFINITY,
        };
    }
    let num_intersects = cell_polys
        .iter()
        .enumerate()
        .map(|(pi, poly)| {
            if pi != ci
                && pi != oci
                && check_other_poly_intersect(&lseg, poly)
            {
                1.0
            } else {
                0.0
            }
        })
        .sum::<f64>();
    VertexPairInfo {
        dist: lseg.len,
        num_intersects,
    }
}

impl CoaGenerator {
    /// Calculates a matrix storing whether two vertices have clear line of sight if in contact range.
    pub fn new(
        cell_polys: &[Poly],
        params: CoaParams,
        phys_contact_generator: &PhysicalContactGenerator,
    ) -> CoaGenerator {
        let num_cells = cell_polys.len();
        let contact_bbs = cell_polys
            .iter()
            .map(|cp| cp.bbox.expand_by(2.0 * params.halfmax_dist))
            .collect::<Vec<BBox>>();
        let contact_matrix = gen_contact_matrix(&contact_bbs);
        let mut dat =
            SymCcVvDat::empty(num_cells, VertexPairInfo::infinity());

        for (ci, poly) in cell_polys.iter().enumerate() {
            for (vi, v) in poly.verts.iter().enumerate() {
                for (ocj, opoly) in
                    cell_polys[(ci + 1)..].iter().enumerate()
                {
                    let oci = (ci + 1) + ocj;
                    for (ovi, ov) in opoly.verts.iter().enumerate() {
                        if !(phys_contact_generator
                            .min_dist_to(ci, vi)
                            < params.too_close_dist_sq
                            || phys_contact_generator
                                .min_dist_to(oci, ovi)
                                < params.too_close_dist_sq)
                            && contact_matrix.get(ci, oci)
                        {
                            let lseg = LineSeg2D::new(v, ov);
                            dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                calc_pair_info(
                                    ci, vi, oci, ovi, lseg,
                                    cell_polys,
                                ),
                            );
                        } else {
                            dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                VertexPairInfo::infinity(),
                            );
                        }
                    }
                }
            }
        }
        let mut r = CoaGenerator {
            dat,
            contact_bbs,
            contact_matrix,
            coa_matrix: vec![],
            cell_polys: cell_polys.to_vec(),
            params,
        };
        r.gen_for_all();
        r
    }

    pub fn update_for(
        &mut self,
        cell_ix: usize,
        cell_poly: Poly,
        phys_contact_generator: &PhysicalContactGenerator,
    ) {
        self.cell_polys[cell_ix] = cell_poly;
        let bb = cell_poly.bbox.expand_by(self.params.halfmax_dist);
        self.contact_bbs[cell_ix] = bb;
        // Update contacts.
        for (oci, obb) in self.contact_bbs.iter().enumerate() {
            if oci != cell_ix {
                self.contact_matrix.set(
                    cell_ix,
                    oci,
                    obb.intersects(&bb),
                )
            }
        }
        for (vi, v) in
            self.cell_polys[cell_ix].verts.iter().enumerate()
        {
            for (ocj, opoly) in
                self.cell_polys[(cell_ix + 1)..].iter().enumerate()
            {
                let oci = (cell_ix + 1) + ocj;
                for (ovi, ov) in opoly.verts.iter().enumerate() {
                    if !(phys_contact_generator
                        .min_dist_to(cell_ix, vi)
                        < self.params.too_close_dist_sq
                        || phys_contact_generator
                            .min_dist_to(oci, ovi)
                            < self.params.too_close_dist_sq)
                        && self.contact_matrix.get(cell_ix, oci)
                    {
                        let lseg = LineSeg2D::new(v, ov);
                        self.dat.set(
                            cell_ix,
                            vi,
                            oci,
                            ovi,
                            calc_pair_info(
                                cell_ix,
                                vi,
                                oci,
                                ovi,
                                lseg,
                                &self.cell_polys,
                            ),
                        );
                    } else {
                        self.dat.set(
                            cell_ix,
                            vi,
                            oci,
                            ovi,
                            VertexPairInfo::infinity(),
                        );
                    }
                }
            }
        }
        self.gen_for(cell_ix);
    }

    fn gen_for_all(&mut self) {
        let num_cells = self.contact_matrix.num_cells;
        for ci in 0..num_cells {
            self.gen_for(ci);
        }
    }

    fn gen_for(&mut self, cell_ix: usize) {
        let num_cells = self.contact_matrix.num_cells;
        let mut x_coas = [0.0f64; NVERTS];
        let CoaParams {
            los_penalty,
            vertex_mag,
            distrib_exp,
            ..
        } = self.params;
        for (vi, x_coa) in x_coas.iter_mut().enumerate() {
            for oci in 0..num_cells {
                if oci != cell_ix {
                    for ovi in 0..NVERTS {
                        let VertexPairInfo {
                            dist,
                            num_intersects,
                        } = self.dat.get(cell_ix, vi, oci, ovi);
                        if num_intersects < f64::INFINITY {
                            let los_factor = 1.0
                                / (num_intersects + 1.0)
                                    .powf(los_penalty);
                            let coa_signal = (distrib_exp * dist)
                                .exp()
                                * los_factor;
                            let additional_signal =
                                vertex_mag * coa_signal;
                            *x_coa += additional_signal;
                        }
                    }
                }
            }
        }
        // let mut max_coa = 0.0;
        // for x_coas in all_x_coas.iter() {
        //     for &x_coa in x_coas.iter() {
        //         if x_coa > max_coa {
        //             max_coa = x_coa;
        //         }
        //     }
        // }
        //println!("{}", max_coa);
        self.coa_matrix[cell_ix] = x_coas;
    }

    pub fn fetch_for(&self, cell_ix: usize) -> [f64; NVERTS] {
        self.coa_matrix[cell_ix]
    }
}
