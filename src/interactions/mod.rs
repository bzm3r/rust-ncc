// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::geometry::{calc_dist_point_to_seg, Bbox};
use crate::math::matrices::{Mat, SymMat};
use crate::math::p2d::P2D;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};

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
pub struct ContactMat {
    dat: SymMat<bool>,
    range: f32,
}

impl ContactMat {
    pub fn get(&self, i: usize, j: usize) -> bool {
        self.dat.get(i, j)
    }

    pub fn num_cells(&self) -> usize {
        self.dat.n
    }
}

/// Stores distance between vertices of cells.
pub struct ContactDistMat(Mat<f32>);

#[allow(unused)]
pub struct CloseEdge {
    cell_ix: usize,
    vert_ix: usize,
    dist: f32,
}

impl ContactDistMat {
    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn close_edges_on_cell(&self, ci: usize, vi: usize, oci: usize) -> Vec<CloseEdge> {
        (0..NVERTS)
            .filter_map(|ovi| {
                let dist = self.0.get(ci * NVERTS + vi, oci * NVERTS + vi);
                if dist < f32::INFINITY {
                    Some(CloseEdge {
                        cell_ix: oci,
                        vert_ix: ovi,
                        dist,
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<CloseEdge>>()
    }

    /// Get edges which contain points close to vertex `vi` on cell `ci`.
    pub fn close_edges(&self, ci: usize, vi: usize) -> Vec<CloseEdge> {
        let mut r = vec![];
        let num_cells = self.0.n_r / NVERTS;
        for oci in 0..num_cells {
            r.append(&mut self.close_edges_on_cell(ci, vi, oci))
        }
        r
    }

    pub fn num_cells(&self) -> usize {
        self.0.n_r / NVERTS
    }
}

/// Calculate a symmetric matrix whose `(i, j)` element is a boolean representing whether cells `i`
/// and `j` are in contact range.
pub fn calc_contact_mat(bboxes: &[Bbox], contact_range: f32) -> ContactMat {
    let num_cells: usize = bboxes.len();
    let bboxes = bboxes
        .iter()
        .map(|bbox| bbox.expand_by(contact_range))
        .collect::<Vec<Bbox>>();
    let mut contact_mat = SymMat::new(num_cells, false);
    for ci in 0..num_cells {
        let this_bbox = &bboxes[ci];
        for oci in (ci + 1)..num_cells {
            contact_mat.set(ci, oci, bboxes[oci].intersects(this_bbox));
        }
    }
    ContactMat {
        dat: contact_mat,
        range: contact_range,
    }
}

/// Calculate distance between vertices of cells to edges of other cells. Rows are cell vertices,
/// and columns are cell edges, defined by the smaller periodic (periodicity given by number of
/// vertices) index of the two vertices which define the edge. Thus, edge `vi` refers to the edge
/// formed by vertices `vi` and `vi + 1 mod NVERTS`.
pub fn calc_contact_dist_mat(
    cell_vcs: &[[P2D; NVERTS]],
    contact_mat: &ContactMat,
) -> ContactDistMat {
    let num_cells = contact_mat.num_cells();
    let mut dist_mat = Mat::new(num_cells * NVERTS, num_cells * NVERTS, f32::INFINITY);
    for ci in 0..num_cells {
        let this_vcs = &cell_vcs[ci];
        for vi in 0..NVERTS {
            for oci in 0..num_cells {
                if contact_mat.get(ci, oci) {
                    let other_vcs = &cell_vcs[oci];
                    for ovi in 0..(NVERTS) {
                        let ovi2 = circ_ix_plus(ovi, NVERTS);
                        let d = calc_dist_point_to_seg(
                            &this_vcs[vi],
                            &other_vcs[ovi],
                            &other_vcs[ovi2],
                        );
                        if d < contact_mat.range || (d - contact_mat.range).abs() < f32::EPSILON {
                            dist_mat.set(ci * NVERTS + vi, oci, d);
                            break;
                        }
                    }
                }
            }
        }
    }
    ContactDistMat(dist_mat)
}

#[derive(Clone, Default, Deserialize, Serialize, Schematize)]
pub struct CloseCellInfo {
    pub num: usize,
    pub indices: [usize; 9],
}

impl CloseCellInfo {
    pub fn insert(&mut self, ix: usize) {
        if self.num < 9 {
            self.indices[self.num] = ix;
            self.num += 1;
        }
    }
}

impl InteractionState {
    pub fn from_contacts(
        contact_dist_mat: &ContactDistMat,
        cil_mat: &CilMat,
    ) -> (Vec<InteractionState>, Vec<CloseCellInfo>) {
        let num_cells = contact_dist_mat.num_cells();
        let mut closest_cells = vec![CloseCellInfo::default(); num_cells];
        let mut interactions = vec![InteractionState::default(); num_cells];
        for ci in 0..num_cells {
            for vi in 0..NVERTS {
                let close_points = contact_dist_mat.close_edges(ci, vi);
                for CloseEdge { cell_ix: oci, .. } in close_points.into_iter() {
                    closest_cells[ci].insert(oci);
                    interactions[ci].x_cils[vi] = cil_mat.get(ci, oci);
                    //let adh_f1 = cells[oci].mech_state.sum_fs[ovi];
                    //let adh_f2 = cells[oci].mech_state.sum_fs[circ_ix_plus(ovi, NVERTS)];
                    interactions[ci].x_adhs[vi] = P2D { x: 0.0, y: 0.0 }; //interactions[ci].x_adhs[vi] + adh_f1 + adh_f2;
                }
            }
        }
        (interactions, closest_cells)
    }
}
