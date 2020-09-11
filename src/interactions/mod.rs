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
pub struct InteractionState {
    pub x_cils: [f32; NVERTS],
    pub x_adhs: [P2D; NVERTS],
    pub x_chemoas: [f32; NVERTS],
    pub x_coas: [f32; NVERTS],
    pub x_bdrys: [f32; NVERTS],
}

#[derive(Clone)]
pub struct CilMat {
    dat: SymMat<f32>,
}

impl CilMat {
    pub fn new(num_cells: usize, default: f32) -> CilMat {
        CilMat {
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
pub struct Contacts {
    dat: SymMat<bool>,
    range: f32,
}

impl Contacts {
    pub fn num_cells(&self) -> usize {
        self.dat.n
    }

    pub fn get_contacts(&self, ci: usize) -> Vec<usize> {
        (0..self.dat.n)
            .filter(|&oci| self.dat.get(ci, oci))
            .collect()
    }

    pub fn in_contact(&self, ci: usize, oci: usize) -> bool {
        self.dat.get(ci, oci)
    }

    /// Calculate distances between vertices of cells in contact.
    pub fn calc_dists(&self, cell_vcs: &[[P2D; NVERTS]]) -> IVDists {
        let num_cells = self.num_cells();
        let mut dist_mat = IVDists::new(num_cells);
        for (ci, vcs) in cell_vcs.iter().enumerate() {
            for (oci, ovcs) in cell_vcs.iter().enumerate() {
                if ci != oci && self.in_contact(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let ovc2 = &ovcs[circ_ix_plus(ovi, NVERTS)];
                            let (t, d) = calc_dist_point_to_seg(vc, ovc, ovc2);
                            match d.partial_cmp(&self.range) {
                                Some(ord) => match ord {
                                    Ordering::Greater => {
                                        continue;
                                    }
                                    _ => {
                                        dist_mat.set(ci, vi, oci, ovi, (t, d));
                                        break;
                                    }
                                },
                                _ => panic!("cannot compare {} and  {}", d, &self.range),
                            }
                        }
                    }
                }
            }
        }
        dist_mat
    }
}

/// Stores inter-vertex distances, if these distances are less than `cutoff`.
#[derive(Clone)]
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

/// Information relevant to calculating CIL and adhesion due to an edge in proximity to the focus
/// vertex.
pub struct CloseEdge {
    /// Cell containing edge close to focus vertex.
    pub cell_ix: usize,
    /// Close edge runs from `vert_ix` to `vert_ix + 1`.
    pub vert_ix: usize,
    /// Let the position of the point on the close edge closest to the focus vertex be denoted `p`,
    /// and the position of the focus vertex be denoted `v`. `delta` is such that `delta + v = p`.
    pub delta: P2D,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub t: f32,
}

impl IVDists {
    pub fn new(num_cells: usize) -> IVDists {
        let stride_v = (num_cells - 1) * NVERTS as usize;
        let stride_c = stride_v * NVERTS as usize;
        IVDists {
            num_cells,
            stride_c,
            stride_v,
            dat: vec![(f32::INFINITY, f32::INFINITY); num_cells * stride_c],
        }
    }

    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn close_edges_on_cell(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        cell_vcs: &Vec<[P2D; 16]>,
        filter: f32,
    ) -> Vec<CloseEdge> {
        let v = cell_vcs[ci][vi];
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
                    Some(CloseEdge {
                        cell_ix: oci,
                        vert_ix: ovi,
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
        cell_vcs: &Vec<[P2D; NVERTS]>,
        filter: f32,
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.num_cells {
            r.append(&mut self.close_edges_on_cell(ci, vi, oci, cell_vcs, filter))
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

    pub fn gen_interactions(
        &self,
        cell_vcs: &Vec<[P2D; NVERTS]>,
        cil_mat: &CilMat,
        adh_const: f32,
        adh_rest: f32,
        adh_criterion: f32,
    ) -> Vec<InteractionState> {
        let num_cells = self.num_cells;
        let mut interactions = vec![InteractionState::default(); num_cells];
        for ci in 0..num_cells {
            for vi in 0..NVERTS {
                for CloseEdge {
                    cell_ix: oci,
                    vert_ix: ovi,
                    delta,
                    t,
                } in self
                    .close_edges(ci, vi, cell_vcs, adh_criterion)
                    .into_iter()
                {
                    interactions[ci].x_cils[vi] = cil_mat.get(ci, oci);
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

/// Determine intercellular contacts within `contact_range`.
pub fn find_possible_contacts(bboxes: &[BBox], contact_range: f32) -> Contacts {
    let num_cells: usize = bboxes.len();
    let bboxes = bboxes
        .iter()
        .map(|bbox| bbox.expand_by(contact_range))
        .collect::<Vec<BBox>>();
    let mut contact_mat = SymMat::new(num_cells, false);
    for (ci, bb) in bboxes.iter().enumerate() {
        for (oxi, obb) in bboxes[(ci + 1)..].iter().enumerate() {
            let intersects = obb.intersects(bb);
            contact_mat.set(ci, ci + 1 + oxi, intersects);
        }
    }
    Contacts {
        dat: contact_mat,
        range: contact_range,
    }
}
