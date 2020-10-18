// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::geometry::{calc_dist_point_to_seg, BBox};
use crate::math::matrices::SymMat;
use crate::math::p2d::P2D;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize)]
pub struct CellInteractions {
    pub x_cals: [f32; NVERTS],
    pub x_cils: [f32; NVERTS],
    pub x_adhs: [P2D; NVERTS],
    pub x_chemoas: [f32; NVERTS],
    pub x_coas: [f32; NVERTS],
    pub x_bdrys: [f32; NVERTS],
}

#[derive(Clone, Default)]
pub struct CrlMat {
    dat: SymMat<f32>,
}

impl CrlMat {
    pub fn new(num_cells: usize, default: f32) -> CrlMat {
        CrlMat {
            dat: SymMat::<f32>::new(num_cells, default),
        }
    }

    pub fn num_cells(&self) -> usize {
        self.dat.n
    }

    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.dat.get(i, j)
    }

    pub fn set(&mut self, i: usize, j: usize, cil: f32) {
        self.dat.set(i, j, cil);
    }
}

/// Stores whether cells are in contact.
#[derive(Clone, Default)]
pub struct InteractionState {
    cell_bbs: Vec<BBox>,
    cell_vcs: Vec<[P2D; NVERTS]>,
    cell_rgtps: Vec<[RgtpState; NVERTS]>,
    contacts: SymMat<bool>,
    iv_dists: IVDists,
    pub cell_interactions: Vec<CellInteractions>,
    contact_range: f32,
    cal_mat: CrlMat,
    cil_mat: CrlMat,
    adh_const: f32,
    adh_criterion: f32,
}

/// Stores intercellular contacts within distance `range`.
impl InteractionState {
    pub fn new(
        cell_bbs: &[BBox],
        cell_vcs: &[[P2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        contact_range: f32,
        cal_mat: CrlMat,
        cil_mat: CrlMat,
        adh_const: f32,
        adh_criterion: f32,
    ) -> InteractionState {
        let num_cells: usize = cell_bbs.len();
        let bboxes = cell_bbs
            .iter()
            .map(|bbox| bbox.expand_by(contact_range))
            .collect::<Vec<BBox>>();
        let mut contacts = SymMat::new(num_cells, false);
        for (ci, bb) in bboxes.iter().enumerate() {
            for (oxi, obb) in bboxes[(ci + 1)..].iter().enumerate() {
                let intersects = obb.intersects(bb);
                contacts.set(ci, ci + 1 + oxi, intersects);
            }
        }
        let iv_dists = IVDists::new(cell_vcs, contact_range, &contacts);
        let cell_interactions = Self::init_cell_interactions(
            cell_vcs,
            cell_rgtps,
            &cal_mat,
            &cil_mat,
            &iv_dists,
            adh_const,
            adh_criterion,
        );
        InteractionState {
            cell_bbs: bboxes.iter().copied().collect(),
            cell_vcs: cell_vcs.iter().copied().collect(),
            cell_rgtps: cell_rgtps.iter().copied().collect(),
            iv_dists,
            contacts,
            contact_range,
            cell_interactions,
            cal_mat,
            cil_mat,
            adh_const,
            adh_criterion,
        }
    }

    pub fn update(&mut self, cell_ix: usize, vcs: &[P2D; NVERTS]) {
        let bbox = BBox::from_points(vcs).expand_by(self.contact_range);
        self.cell_bbs[cell_ix] = bbox;
        self.cell_vcs[cell_ix] = *vcs;
        for (oci, obb) in self.cell_bbs.iter().enumerate() {
            if oci != cell_ix {
                let intersects = obb.intersects(&bbox);
                self.contacts.set(cell_ix, oci, intersects);
            }
        }
        self.update_dists(cell_ix);
    }

    fn update_dists(&mut self, cell_ix: usize) {
        let vcs = self.cell_vcs[cell_ix];
        for (oci, ovcs) in self.cell_vcs.iter().enumerate() {
            if cell_ix != oci && self.contacts.get(cell_ix, oci) {
                for (vi, vc) in vcs.iter().enumerate() {
                    for (ovi, ovc) in ovcs.iter().enumerate() {
                        let ovc2 = &ovcs[circ_ix_plus(ovi, NVERTS)];
                        let (t, d) = calc_dist_point_to_seg(vc, ovc, ovc2);
                        match d.partial_cmp(&self.contact_range) {
                            Some(ord) => match ord {
                                Ordering::Greater => {
                                    continue;
                                }
                                _ => {
                                    self.iv_dists.set(cell_ix, vi, oci, ovi, (t, d));
                                    break;
                                }
                            },
                            _ => panic!("cannot compare {} and  {}", d, &self.contact_range),
                        }
                    }
                }
            }
        }
    }

    pub fn update_cell_interactions(&mut self) {
        self.cell_interactions = Self::init_cell_interactions(
            &self.cell_vcs,
            &self.cell_rgtps,
            &self.cal_mat,
            &self.cil_mat,
            &self.iv_dists,
            self.adh_const,
            self.adh_criterion,
        );
    }

    pub fn get_contacts(&self, ci: usize) -> Vec<usize> {
        (0..self.cell_bbs.len())
            .filter(|&oci| self.contacts.get(ci, oci))
            .collect()
    }

    pub fn get_contact_polys(&self, ci: usize) -> Vec<(BBox, [P2D; NVERTS])> {
        (0..self.cell_bbs.len())
            .filter(|&oci| self.contacts.get(ci, oci))
            .map(|oci| {
                let vcs = self.cell_vcs[oci];
                let bb = BBox::from_points(&vcs);
                (bb, vcs)
            })
            .collect()
    }

    pub fn init_cell_interactions(
        cell_vcs: &[[P2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        cal_mat: &CrlMat,
        cil_mat: &CrlMat,
        iv_dists: &IVDists,
        adh_const: f32,
        adh_criterion: f32,
    ) -> Vec<CellInteractions> {
        let num_cells = cell_vcs.len();
        let mut interactions = vec![CellInteractions::default(); num_cells];
        for ci in 0..num_cells {
            for vi in 0..NVERTS {
                for CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    crl,
                    delta,
                    t,
                } in iv_dists
                    .close_edges(ci, vi, cell_vcs, cell_rgtps, adh_criterion)
                    .into_iter()
                {
                    match crl {
                        CrlEffect::Cal => interactions[ci].x_cals[vi] = cal_mat.get(ci, oci),
                        CrlEffect::Cil => interactions[ci].x_cils[vi] = cil_mat.get(ci, oci),
                    };

                    let d = delta.mag();
                    let adh = if d.abs() < 1e-3 {
                        P2D::default()
                    } else {
                        adh_const * d * delta.unitize()
                    };
                    interactions[oci].x_adhs[ovi] =
                        interactions[oci].x_adhs[ovi] + -1.0 * (1.0 - t) * adh;
                    interactions[oci].x_adhs[circ_ix_plus(ovi, NVERTS)] =
                        interactions[oci].x_adhs[circ_ix_plus(ovi, NVERTS)] + -1.0 * t * adh;
                    interactions[ci].x_adhs[vi] = interactions[ci].x_adhs[vi] + adh;
                }
            }
        }
        interactions
    }
}

/// Stores inter-vertex distances, if these distances are less than `cutoff`.
#[derive(Default, Clone)]
pub struct IVDists {
    num_cells: usize,
    /// Number of elements between data of different cells.
    stride_c: usize,
    /// Number of elements between data of different vertices of a cell.
    stride_v: usize,
    /// Stores `(t, d)`, where `(t, d)` is information returned by
    /// `crate::math::geometry::calc_dist_point_to_seg`.
    dat: Vec<(f32, f32)>,
}

pub type RgtpState = f32;

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
    pub delta: P2D,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub t: f32,
}

impl IVDists {
    pub fn empty(num_cells: usize) -> IVDists {
        let stride_v = (num_cells - 1) * NVERTS as usize;
        let stride_c = stride_v * NVERTS as usize;
        IVDists {
            num_cells,
            stride_c,
            stride_v,
            dat: vec![(f32::INFINITY, f32::INFINITY); num_cells * stride_c],
        }
    }

    /// Calculate distances between vertices of cells in contact.
    pub fn new(cell_vcs: &[[P2D; NVERTS]], contact_range: f32, contacts: &SymMat<bool>) -> IVDists {
        let num_cells = cell_vcs.len();
        let mut dist_mat = IVDists::empty(num_cells);
        for (ci, vcs) in cell_vcs.iter().enumerate() {
            for (oci, ovcs) in cell_vcs.iter().enumerate() {
                if ci != oci && contacts.get(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let ovc2 = &ovcs[circ_ix_plus(ovi, NVERTS)];
                            let (t, d) = calc_dist_point_to_seg(vc, ovc, ovc2);
                            match d.partial_cmp(&contact_range) {
                                Some(ord) => match ord {
                                    Ordering::Greater => {
                                        continue;
                                    }
                                    _ => {
                                        dist_mat.set(ci, vi, oci, ovi, (t, d));
                                        break;
                                    }
                                },
                                _ => panic!("cannot compare {} and  {}", d, &contact_range),
                            }
                        }
                    }
                }
            }
        }
        dist_mat
    }

    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn close_edges_on_cell(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        cell_vcs: &[[P2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        filter: f32,
    ) -> Vec<CloseEdge> {
        let v = cell_vcs[ci][vi];
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| {
                #[cfg(debug_assertions)]
                self.check_indices(ci, vi, oci, ovi);
                let (t, d) = self.get(ci, vi, oci, ovi);
                if d < filter {
                    let p0 = cell_vcs[oci][ovi];
                    let p1 = cell_vcs[oci][circ_ix_plus(ovi, NVERTS)];
                    let p = t * (p1 - p0) + p0;
                    let delta = p - v;

                    let edge_rgtp =
                        (cell_rgtps[oci][ovi] + cell_rgtps[oci][circ_ix_plus(ovi, NVERTS)]) / 2.0;
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
        cell_vcs: &[[P2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        filter: f32,
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.num_cells {
            r.append(&mut self.close_edges_on_cell(ci, vi, oci, cell_vcs, cell_rgtps, filter))
        }
        r
    }

    #[cfg(debug_assertions)]
    pub fn check_indices(&self, ci: usize, vi: usize, oci: usize, ovi: usize) {
        if ci > self.num_cells - 1 {
            panic!("{} cells tracked, received ci: {}", self.num_cells, ci);
        }

        if vi > NVERTS as usize {
            panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
        }

        if oci > self.num_cells {
            panic!("{} cells tracked, received oci: {}", self.num_cells, oci);
        }

        if ovi > NVERTS as usize {
            panic!("{} vertices tracked, received vi: {}", NVERTS, vi);
        }
    }

    pub fn calc_ix(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> Option<usize> {
        #[cfg(debug_assertions)]
        self.check_indices(ci, vi, oci, ovi);
        match oci.cmp(&ci) {
            Ordering::Greater => {
                Some(ci * self.stride_c + vi * self.stride_v + (oci - 1) * (NVERTS as usize) + ovi)
            }
            Ordering::Less => {
                Some(ci * self.stride_c + vi * self.stride_v + oci * (NVERTS as usize) + ovi)
            }
            Ordering::Equal => None,
        }
    }

    pub fn set(&mut self, ci: usize, vi: usize, oci: usize, ovi: usize, td: (f32, f32)) {
        if let Some(i) = self.calc_ix(ci, vi, oci, ovi) {
            self.dat[i] = td;
        }
    }

    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> (f32, f32) {
        if let Some(i) = self.calc_ix(ci, vi, oci, ovi) {
            self.dat[i]
        } else {
            (f32::INFINITY, f32::INFINITY)
        }
    }
}
