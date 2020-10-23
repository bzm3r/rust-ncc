// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::geometry::{
    calc_dist_point_to_seg, ls_intersects_poly, ls_self_intersects_poly, BBox, LineSeg,
};
use crate::math::matrices::{CvCvDat, SymCcDat, SymCcVvDat};
use crate::math::p2d::V2D;
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize)]
pub struct CellInteractions {
    pub x_cals: [f32; NVERTS],
    pub x_cils: [f32; NVERTS],
    pub x_adhs: [V2D; NVERTS],
    pub x_chemoas: [f32; NVERTS],
    pub x_coas: [f32; NVERTS],
    pub x_bdrys: [f32; NVERTS],
}

#[derive(Clone, Default)]
pub struct ContactInfo {}

/// Stores whether cells are in contact.
#[derive(Clone)]
pub struct InteractionState {
    cell_vcs: Vec<[V2D; NVERTS]>,
    cell_rgtps: Vec<[RgtpState; NVERTS]>,
    cal_mat: SymCcDat<f32>,
    cil_mat: SymCcDat<f32>,
    coa_range: f32,
    crl_range: f32,
    coa_bbs: Vec<BBox>,
    crl_bbs: Vec<BBox>,
    coa_contacts: SymCcDat<bool>,
    crl_contacts: SymCcDat<bool>,
    crl_dists: CrlDists,
    coa_dists: CoaDists,
    pub cell_interactions: Vec<CellInteractions>,
    adh_const: f32,
    adh_criterion: f32,
}

/// Stores intercellular contacts within distance `range`.
impl InteractionState {
    pub fn new(
        cell_bbs: &[BBox],
        cell_vcs: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        crl_range: f32,
        coa_range: f32,
        cal_mat: SymCcDat<f32>,
        cil_mat: SymCcDat<f32>,
        adh_const: f32,
        adh_criterion: f32,
    ) -> InteractionState {
        let coa_bbs = cell_bbs
            .iter()
            .map(|bb| bb.expand_by(coa_range))
            .collect::<Vec<BBox>>();
        let crl_bbs = cell_bbs
            .iter()
            .map(|bb| bb.expand_by(coa_range))
            .collect::<Vec<BBox>>();
        let num_cells = cell_vcs.len();
        let mut coa_contacts = SymCcDat::new(num_cells, false);
        let mut crl_contacts = SymCcDat::new(num_cells, false);
        for (ci, (coa_bb, crl_bb)) in coa_bbs.iter().zip(crl_bbs.iter()).enumerate() {
            for (oxi, (ocoa, ocrl)) in coa_bbs[(ci + 1)..]
                .iter()
                .zip(crl_bbs[(ci + 1)..].iter())
                .enumerate()
            {
                coa_contacts.set(ci, ci + 1 + oxi, ocoa.intersects(coa_bb));
                crl_contacts.set(ci, ci + 1 + oxi, ocrl.intersects(crl_bb));
            }
        }
        let coa_dists = CoaDists::new(cell_vcs, &coa_contacts);
        let crl_dists = CrlDists::new(cell_vcs, crl_range, &crl_contacts);
        let cell_interactions = Self::init_cell_interactions(
            cell_vcs,
            cell_rgtps,
            &cal_mat,
            &cil_mat,
            &coa_dists,
            &crl_dists,
            adh_const,
            adh_criterion,
        );
        InteractionState {
            cell_vcs: cell_vcs.iter().copied().collect(),
            cell_rgtps: cell_rgtps.iter().copied().collect(),
            crl_dists,
            coa_dists,
            cell_interactions,
            cal_mat,
            cil_mat,
            coa_range,
            crl_range,
            coa_bbs,
            crl_bbs,
            coa_contacts,
            adh_const,
            adh_criterion,
            crl_contacts,
        }
    }

    pub fn update_contacts(&mut self, cell_ix: usize, vcs: &[V2D; NVERTS]) {
        let bb = BBox::from_points(vcs);
        let coa_bb = bb.expand_by(self.coa_range);
        let crl_bb = bb.expand_by(self.crl_range);
        self.coa_bbs[cell_ix] = coa_bb;
        self.crl_bbs[cell_ix] = crl_bb;
        for (oci, (ocoa, ocrl)) in self.coa_bbs.iter().zip(self.crl_bbs.iter()).enumerate() {
            if oci != cell_ix {
                self.coa_contacts
                    .set(cell_ix, oci, ocoa.intersects(&coa_bb));
                self.crl_contacts
                    .set(cell_ix, oci, ocrl.intersects(&coa_bb));
            }
        }
    }

    pub fn get_coa_contacts(&self, ci: usize) -> Vec<usize> {
        (0..self.coa_bbs.len())
            .filter(|&oci| self.coa_contacts.get(ci, oci))
            .collect()
    }

    pub fn get_crl_contacts(&self, ci: usize) -> Vec<usize> {
        (0..self.crl_bbs.len())
            .filter(|&oci| self.crl_contacts.get(ci, oci))
            .collect()
    }

    pub fn get_coa_contact_bbs(&self, ci: usize) -> Vec<BBox> {
        self.coa_bbs
            .iter()
            .enumerate()
            .filter_map(|(oci, &obb)| {
                if self.coa_contacts.get(ci, oci) {
                    Some(obb)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_crl_contact_bbs(&self, ci: usize) -> Vec<BBox> {
        self.crl_bbs
            .iter()
            .enumerate()
            .filter_map(|(oci, &obb)| {
                if self.crl_contacts.get(ci, oci) {
                    Some(obb)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn update(&mut self, cell_ix: usize, vcs: &[V2D; NVERTS]) {
        self.update_contacts(cell_ix, vcs);
        self.update_dists(cell_ix);
    }

    fn update_dists(&mut self, ci: usize) {
        self.crl_dists
            .update(ci, &self.cell_vcs, self.crl_range, &self.crl_contacts);
        self.coa_dists
            .update(ci, &self.cell_vcs, &self.coa_contacts);
    }

    pub fn update_cell_interactions(&mut self) {
        self.cell_interactions = Self::init_cell_interactions(
            &self.cell_vcs,
            &self.cell_rgtps,
            &self.cal_mat,
            &self.cil_mat,
            &self.coa_dists,
            &self.crl_dists,
            self.adh_const,
            self.adh_criterion,
        );
    }

    pub fn init_cell_interactions(
        cell_vcs: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        cal_mat: &SymCcDat<f32>,
        cil_mat: &SymCcDat<f32>,
        coa_dists: &CoaDists,
        crl_dists: &CrlDists,
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
                } in crl_dists
                    .close_edges(ci, vi, cell_vcs, cell_rgtps, adh_criterion)
                    .into_iter()
                {
                    match crl {
                        CrlEffect::Cal => interactions[ci].x_cals[vi] = cal_mat.get(ci, oci),
                        CrlEffect::Cil => interactions[ci].x_cils[vi] = cil_mat.get(ci, oci),
                    };

                    let d = delta.mag();
                    let adh = if d.abs() < 1e-3 {
                        V2D::default()
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
        for ci in 0..num_cells {
            for oci in (ci + 1)..num_cells {
                for vi in 0..NVERTS {
                    for ovi in 0..NVERTS {
                        let (d, los) = coa_dists.get(ci, vi, oci, ovi);
                        if los != f32::INFINITY && d != f32::INFINITY {}
                    }
                }
            }
        }
        interactions
    }
}

#[derive(Clone)]
pub struct CrlDists {
    dat: CvCvDat<(f32, f32)>,
}

impl CrlDists {
    /// Calculate distances between vertices of cells in contact.
    pub fn new(cell_vcs: &[[V2D; NVERTS]], crl_range: f32, contacts: &SymCcDat<bool>) -> CrlDists {
        let num_cells = cell_vcs.len();
        let mut dat = CvCvDat::empty(num_cells, (f32::INFINITY, f32::INFINITY));
        for (ci, vcs) in cell_vcs.iter().enumerate() {
            for (oci, ovcs) in cell_vcs.iter().enumerate() {
                if ci != oci && contacts.get(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let ovc2 = &ovcs[circ_ix_plus(ovi, NVERTS)];
                            let (t, d) = calc_dist_point_to_seg(vc, ovc, ovc2);
                            match d.partial_cmp(&crl_range) {
                                Some(ord) => match ord {
                                    Ordering::Greater => {
                                        continue;
                                    }
                                    _ => {
                                        dat.set(ci, vi, oci, ovi, (t, d));
                                        break;
                                    }
                                },
                                _ => panic!("cannot compare {} and  {}", d, &crl_range),
                            }
                        }
                    }
                }
            }
        }
        CrlDists { dat }
    }

    /// Get edges containing points on cell `oci` which are close to vertex `vi` on cell `ci`.
    pub fn close_edges_on_cell(
        &self,
        ci: usize,
        vi: usize,
        oci: usize,
        cell_vcs: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        filter: f32,
    ) -> Vec<CloseEdge> {
        let v = cell_vcs[ci][vi];
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| {
                let (t, d) = self.dat.get(ci, vi, oci, ovi);
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
        cell_vcs: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        filter: f32,
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.dat.num_cells {
            r.append(&mut self.close_edges_on_cell(ci, vi, oci, cell_vcs, cell_rgtps, filter))
        }
        r
    }

    pub fn update(
        &mut self,
        ci: usize,
        cell_vcs: &[[V2D; NVERTS]],
        crl_range: f32,
        crl_contacts: &SymCcDat<bool>,
    ) {
        let vcs = &cell_vcs[ci];
        for (oci, ovcs) in cell_vcs.iter().enumerate() {
            if ci != oci && crl_contacts.get(ci, oci) {
                for (vi, vc) in vcs.iter().enumerate() {
                    for (ovi, ovc) in ovcs.iter().enumerate() {
                        let ovc2 = &ovcs[circ_ix_plus(ovi, NVERTS)];
                        let (t, d) = calc_dist_point_to_seg(vc, ovc, ovc2);
                        match d.partial_cmp(&crl_range) {
                            Some(ord) => match ord {
                                Ordering::Greater => {
                                    continue;
                                }
                                _ => {
                                    self.dat.set(ci, vi, oci, ovi, (t, d));
                                    break;
                                }
                            },
                            _ => panic!("cannot compare {} and  {}", d, crl_range),
                        }
                    }
                }
            }
        }
    }
}

#[derive(Clone)]
pub struct CoaDists {
    dat: SymCcVvDat<(f32, f32)>,
}

impl CoaDists {
    /// Calculates a matrix storing whether two vertices have clear line of sight if in contact range.
    pub fn new(cell_vcs: &[[V2D; NVERTS]], coa_contacts: &SymCcDat<bool>) -> CoaDists {
        let num_cells = cell_vcs.len();
        let mut dat = SymCcVvDat::empty(num_cells, (f32::INFINITY, f32::INFINITY));
        for (ci, vcs) in cell_vcs.iter().enumerate() {
            for (ocj, ovcs) in cell_vcs[(ci + 1)..].iter().enumerate() {
                let oci = ci + ocj;
                if ci != oci && coa_contacts.get(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let lseg = LineSeg::new(vc, ovc);
                            if ls_self_intersects_poly(vi, vcs, &lseg)
                                || ls_self_intersects_poly(ovi, ovcs, &lseg)
                            {
                                dat.set(ci, vi, oci, ovi, (f32::INFINITY, f32::INFINITY));
                            } else {
                                dat.set(
                                    ci,
                                    vi,
                                    oci,
                                    ovi,
                                    (
                                        (vc - ovc).mag(),
                                        cell_vcs
                                            .iter()
                                            .enumerate()
                                            .map(|(pi, poly)| {
                                                if pi != ci
                                                    && pi != oci
                                                    && ls_intersects_poly(&lseg, poly)
                                                {
                                                    1.0
                                                } else {
                                                    0.0
                                                }
                                            })
                                            .sum::<f32>(),
                                    ),
                                )
                            }
                        }
                    }
                }
            }
        }
        CoaDists { dat }
    }

    pub fn update(&mut self, ci: usize, cell_vcs: &[[V2D; NVERTS]], coa_contacts: &SymCcDat<bool>) {
        let vcs = &cell_vcs[ci];
        for (ocj, ovcs) in cell_vcs[(ci + 1)..].iter().enumerate() {
            let oci = ci + ocj;
            if ci != oci && coa_contacts.get(ci, oci) {
                for (vi, vc) in vcs.iter().enumerate() {
                    for (ovi, ovc) in ovcs.iter().enumerate() {
                        let lseg = LineSeg::new(vc, ovc);
                        if ls_self_intersects_poly(vi, vcs, &lseg)
                            || ls_self_intersects_poly(ovi, ovcs, &lseg)
                        {
                            self.dat
                                .set(ci, vi, oci, ovi, (f32::INFINITY, f32::INFINITY));
                        } else {
                            self.dat.set(
                                ci,
                                vi,
                                oci,
                                ovi,
                                (
                                    (vc - ovc).mag(),
                                    cell_vcs
                                        .iter()
                                        .enumerate()
                                        .map(|(pi, poly)| {
                                            if pi != ci
                                                && pi != oci
                                                && ls_intersects_poly(&lseg, poly)
                                            {
                                                1.0
                                            } else {
                                                0.0
                                            }
                                        })
                                        .sum::<f32>(),
                                ),
                            )
                        }
                    }
                }
            }
        }
    }

    pub fn get(&self, ci: usize, vi: usize, oci: usize, ovi: usize) -> (f32, f32) {
        self.dat.get(ci, vi, oci, ovi)
    }
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
    pub delta: V2D,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub t: f32,
}
