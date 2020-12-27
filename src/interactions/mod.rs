// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

pub mod dat2d;
pub mod dat4d;
pub mod symdat2d;
pub mod symdat4d;
mod utils;

use crate::interactions::dat4d::CvCvDat;
use crate::interactions::symdat2d::SymCcDat;
use crate::interactions::symdat4d::SymCcVvDat;
use crate::math::geometry::{
    calc_dist_point_to_seg, is_point_in_poly,
    is_point_in_poly_no_bb_check, ls_intersects_poly,
    ls_self_intersects_poly, BBox, LineSeg,
};
use crate::math::v2d::V2d;
use crate::parameters::{
    BdryParams, ChemAttrParams, CoaParams, InteractionParams,
    PhysicalContactParams,
};
use crate::utils::circ_ix_plus;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(
    Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize,
)]
pub struct CellInteractions {
    pub x_cals: [f32; NVERTS],
    pub x_cils: [f32; NVERTS],
    pub x_adhs: [V2d; NVERTS],
    pub x_chem_attrs: [f32; NVERTS],
    pub x_coas: [f32; NVERTS],
    pub x_bdrys: [f32; NVERTS],
}

#[derive(Clone, Default)]
pub struct ContactInfo {}

/// Generates interaction related factors.
#[derive(Clone)]
pub struct InteractionGenerator {
    /// Vertex coordinates, per cell, for all cells in the simulation.
    cell_vcs: Vec<[V2d; NVERTS]>,
    cell_rgtps: Vec<[RgtpState; NVERTS]>,
    /// Generates CIL/CAL related interaction information. In other words,
    /// interactions that require cells to engage in physical contact.
    phys_contact_generator: PhysicalContactGenerator,
    coa_generator: CoaGenerator,
    chem_attr_generator: ChemAttrGenerator,
    bdry_generator: BdryGenerator,
}

/// Stores intercellular contacts within distance `range`.
impl InteractionGenerator {
    pub fn new(
        cell_vcs: &[[V2d; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        parameters: InteractionParams,
    ) -> InteractionGenerator {
        let cell_bbs = cell_vcs
            .iter()
            .map(|vcs| BBox::from_points(vcs))
            .collect::<Vec<BBox>>();
        let InteractionParams {
            phys_contact: phys_contact_params,
            coa: coa_params,
            chem_attr: chem_attr_params,
            bdry: bdry_params,
        } = parameters;
        let coa_generator =
            CoaGenerator::new(&cell_bbs, cell_vcs, coa_params);
        let phys_contact_generator = PhysicalContactGenerator::new(
            &cell_bbs,
            cell_vcs,
            phys_contact_params,
        );
        let chem_attr_generator =
            ChemAttrGenerator::new(chem_attr_params);
        let bdry_generator = BdryGenerator::new(bdry_params);
        InteractionGenerator {
            cell_vcs: cell_vcs.iter().copied().collect(),
            cell_rgtps: cell_rgtps.iter().copied().collect(),
            phys_contact_generator,
            coa_generator,
            chem_attr_generator,
            bdry_generator,
        }
    }

    pub fn update(&mut self, cell_ix: usize, vcs: &[V2d; NVERTS]) {
        self.cell_vcs[cell_ix] = *vcs;
        let bb = BBox::from_points(vcs);
        self.coa_generator.update(cell_ix, &bb, vcs, &self.cell_vcs);
        self.phys_contact_generator.update(
            cell_ix,
            &bb,
            vcs,
            &self.cell_vcs,
        );
    }

    pub fn generate(&self) -> Vec<CellInteractions> {
        let num_cells = self.cell_vcs.len();
        let (r_adhs, r_cals, r_cils) = self
            .phys_contact_generator
            .generate(&self.cell_vcs, &self.cell_rgtps);
        let r_coas = self.coa_generator.generate();
        let r_chemoas =
            self.chem_attr_generator.generate(&self.cell_vcs);
        let r_bdrys = self.bdry_generator.generate(&self.cell_vcs);
        (0..num_cells)
            .map(|ci| CellInteractions {
                x_cals: r_cals[ci],
                x_cils: r_cils[ci],
                x_adhs: r_adhs[ci],
                x_chem_attrs: r_chemoas[ci],
                x_coas: r_coas[ci],
                x_bdrys: r_bdrys[ci],
            })
            .collect()
    }

    pub fn get_physical_contacts(&self, ci: usize) -> Vec<usize> {
        let num_cells = self.cell_vcs.len();
        (0..num_cells)
            .filter(|&oci| {
                self.phys_contact_generator.contacts.get(ci, oci)
            })
            .collect()
    }

    pub fn get_physical_contact_polys(
        &self,
        ci: usize,
    ) -> Vec<(BBox, [V2d; NVERTS])> {
        let contacts = self.get_physical_contacts(ci);
        contacts
            .into_iter()
            .map(|oci| {
                (
                    self.phys_contact_generator.contact_bbs[oci],
                    self.cell_vcs[oci],
                )
            })
            .collect()
    }
}

pub type ClosestDistLoc = f32;
pub type ClosestDist = f32;

/// Generates CIL/CAL related interaction information. In other words,
/// interactions that require cells to engage in physical contact.
#[derive(Clone)]
pub struct PhysicalContactGenerator {
    ///
    dat: CvCvDat<(ClosestDistLoc, ClosestDist)>,
    contact_bbs: Vec<BBox>,
    contacts: SymCcDat<bool>,
    params: Option<PhysicalContactParams>,
}

impl PhysicalContactGenerator {
    /// Calculate distances between vertices of cells in contact.
    pub fn new(
        cell_bbs: &[BBox],
        cell_vcs: &[[V2d; NVERTS]],
        params: Option<PhysicalContactParams>,
    ) -> PhysicalContactGenerator {
        let num_cells = cell_bbs.len();
        let mut dat =
            CvCvDat::empty(num_cells, (f32::INFINITY, f32::INFINITY));
        let mut contact_bbs = vec![];
        let mut contacts = SymCcDat::new(num_cells, false);
        if let Some(params) = &params {
            cell_bbs.iter().for_each(|bb| {
                contact_bbs.push(bb.expand_by(params.range))
            });
            for (ci, bb) in contact_bbs.iter().enumerate() {
                for (oxi, obb) in
                    contact_bbs[(ci + 1)..].iter().enumerate()
                {
                    contacts.set(
                        ci,
                        ci + 1 + oxi,
                        obb.intersects(bb),
                    );
                }
            }
            for (ci, vcs) in cell_vcs.iter().enumerate() {
                for (oci, ovcs) in cell_vcs.iter().enumerate() {
                    if ci != oci && contacts.get(ci, oci) {
                        for (vi, vc) in vcs.iter().enumerate() {
                            for (ovi, ovc) in ovcs.iter().enumerate()
                            {
                                let ovc2 =
                                    &ovcs[circ_ix_plus(ovi, NVERTS)];
                                let (t, d) = calc_dist_point_to_seg(
                                    vc, ovc, ovc2,
                                );
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
        cell_vcs: &[[V2d; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let v = cell_vcs[ci][vi];
        let v_rgtp = cell_rgtps[ci][vi];
        (0..NVERTS)
            .filter_map(|ovi| {
                let (t, d) = self.dat.get(ci, vi, oci, ovi);
                if d < f32::INFINITY {
                    let p0 = cell_vcs[oci][ovi];
                    let p1 = cell_vcs[oci][circ_ix_plus(ovi, NVERTS)];
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
        cell_vcs: &[[V2d; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
    ) -> Vec<CloseEdge> {
        let mut r = vec![];
        for oci in 0..self.dat.num_cells {
            r.append(&mut self.close_edges_on_cell(
                ci, vi, oci, cell_vcs, cell_rgtps,
            ))
        }
        r
    }

    pub fn update(
        &mut self,
        ci: usize,
        bb: &BBox,
        vcs: &[V2d; NVERTS],
        cell_vcs: &[[V2d; NVERTS]],
    ) {
        if let Some(params) = &self.params {
            let bb = bb.expand_by(params.range);
            self.contact_bbs[ci] = bb;
            for (oxi, obb) in
                self.contact_bbs[(ci + 1)..].iter().enumerate()
            {
                self.contacts.set(
                    ci,
                    ci + 1 + oxi,
                    obb.intersects(&bb),
                );
            }
            for (oci, ovcs) in cell_vcs.iter().enumerate() {
                if ci != oci && self.contacts.get(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let ovc2 =
                                &ovcs[circ_ix_plus(ovi, NVERTS)];
                            let (t, d) =
                                calc_dist_point_to_seg(vc, ovc, ovc2);
                            match d.partial_cmp(&params.range) {
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
                                    d, &params.range
                                ),
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn generate(
        &self,
        cell_vcs: &[[V2d; NVERTS]],
        cell_rgtps: &[[f32; NVERTS]],
    ) -> (Vec<[V2d; NVERTS]>, Vec<[f32; NVERTS]>, Vec<[f32; NVERTS]>)
    {
        let num_cells = self.contacts.num_cells;
        let mut r_adhs = vec![[V2d::default(); NVERTS]; num_cells];
        let mut r_cals = vec![[0.0f32; NVERTS]; num_cells];
        let mut r_cils = vec![[0.0f32; NVERTS]; num_cells];
        if let Some(p) = &self.params {
            for ci in 0..num_cells {
                let x_cals = &mut r_cals[ci];
                let x_cils = &mut r_cils[ci];
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
                                    ci, vi, cell_vcs, cell_rgtps,
                                )
                                .into_iter()
                            {
                                match crl {
                                    CrlEffect::Cal => {
                                        x_cals[vi] = p.cal_mag
                                    }
                                    CrlEffect::Cil => {
                                        x_cils[vi] = p.cil_mag
                                    }
                                };

                                let d = delta.mag();
                                let adh = if d.abs() < 1e-3 {
                                    V2d::default()
                                } else {
                                    p.adh_mag * d * delta.unitize()
                                };
                                r_adhs[oci][ovi] = r_adhs[oci][ovi]
                                    + -1.0 * (1.0 - t) * adh;
                                r_adhs[oci]
                                    [circ_ix_plus(ovi, NVERTS)] =
                                    r_adhs[oci]
                                        [circ_ix_plus(ovi, NVERTS)]
                                        + -1.0 * t * adh;
                                r_adhs[ci][vi] = r_adhs[ci][vi] + adh;
                            }
                        }
                    }
                }
            }
        }
        (r_adhs, r_cals, r_cils)
    }
}

#[derive(Clone)]
pub struct CoaGenerator {
    dat: SymCcVvDat<(f32, f32)>,
    contact_bbs: Vec<BBox>,
    contacts: SymCcDat<bool>,
    params: Option<CoaParams>,
}

impl CoaGenerator {
    /// Calculates a matrix storing whether two vertices have clear line of sight if in contact range.
    pub fn new(
        cell_bbs: &[BBox],
        cell_vcs: &[[V2d; NVERTS]],
        params: Option<CoaParams>,
    ) -> CoaGenerator {
        let num_cells = cell_bbs.len();
        let mut contact_bbs = vec![];
        let mut contacts = SymCcDat::new(num_cells, false);
        let mut dat = SymCcVvDat::empty(
            num_cells,
            (f32::INFINITY, f32::INFINITY),
        );
        if let Some(params) = &params {
            cell_bbs.iter().for_each(|bb| {
                contact_bbs.push(bb.expand_by(params.range))
            });
            for (ci, bb) in contact_bbs.iter().enumerate() {
                for (oxi, obb) in
                    contact_bbs[(ci + 1)..].iter().enumerate()
                {
                    contacts.set(
                        ci,
                        ci + 1 + oxi,
                        obb.intersects(bb),
                    );
                }
            }
            for (ci, vcs) in cell_vcs.iter().enumerate() {
                for (ocj, ovcs) in
                    cell_vcs[(ci + 1)..].iter().enumerate()
                {
                    let oci = ci + ocj;
                    if ci != oci && contacts.get(ci, oci) {
                        for (vi, vc) in vcs.iter().enumerate() {
                            for (ovi, ovc) in ovcs.iter().enumerate()
                            {
                                let lseg = LineSeg::new(vc, ovc);
                                if ls_self_intersects_poly(
                                    vi, vcs, &lseg,
                                ) || ls_self_intersects_poly(
                                    ovi, ovcs, &lseg,
                                ) {
                                    dat.set(
                                        ci,
                                        vi,
                                        oci,
                                        ovi,
                                        (
                                            f32::INFINITY,
                                            f32::INFINITY,
                                        ),
                                    );
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
        vcs: &[V2d; NVERTS],
        cell_vcs: &[[V2d; NVERTS]],
    ) {
        if let Some(params) = &self.params {
            let bb = bb.expand_by(params.range);
            self.contact_bbs[ci] = bb;
            for (oci, obb) in self.contact_bbs.iter().enumerate() {
                if oci != ci {
                    self.contacts.set(ci, oci, obb.intersects(&bb))
                }
            }
            for (ocj, ovcs) in cell_vcs[(ci + 1)..].iter().enumerate()
            {
                let oci = ci + ocj;
                if ci != oci && self.contacts.get(ci, oci) {
                    for (vi, vc) in vcs.iter().enumerate() {
                        for (ovi, ovc) in ovcs.iter().enumerate() {
                            let lseg = LineSeg::new(vc, ovc);
                            if ls_self_intersects_poly(vi, vcs, &lseg)
                                || ls_self_intersects_poly(
                                    ovi, ovcs, &lseg,
                                )
                            {
                                self.dat.set(
                                    ci,
                                    vi,
                                    oci,
                                    ovi,
                                    (f32::INFINITY, f32::INFINITY),
                                );
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
    }

    pub fn generate(&self) -> Vec<[f32; NVERTS]> {
        let num_cells = self.contacts.num_cells;
        let mut r = vec![[0.0f32; NVERTS]; num_cells];
        if let Some(p) = &self.params {
            for (ci, x_coas) in r.iter_mut().enumerate() {
                for (vi, x_coa) in x_coas.iter_mut().enumerate() {
                    for oci in 0..num_cells {
                        if oci != ci {
                            for ovi in 0..NVERTS {
                                let (num_intersects, dist) =
                                    self.dat.get(ci, vi, oci, ovi);
                                *x_coa = p.mag
                                    * (p.distrib_exp * dist).exp()
                                    / (num_intersects + 1.0)
                                        .powf(p.los_penalty);
                            }
                        }
                    }
                }
            }
        }
        r
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
    pub delta: V2d,
    /// Let the position of `vert_ix` be `p0`, and the position of `vert_ix + 1` be `p1`. Let `p`
    /// be the point on the close edge closest to the focus vertex. Then, `t` is such that
    /// `(p1 - p0)*t + p0 = p`.
    pub t: f32,
}

#[derive(Clone)]
pub struct ChemAttrGenerator {
    params: Option<ChemAttrParams>,
}

impl ChemAttrGenerator {
    pub fn new(params: Option<ChemAttrParams>) -> ChemAttrGenerator {
        ChemAttrGenerator { params }
    }

    pub fn generate(
        &self,
        cell_vcs: &[[V2d; NVERTS]],
    ) -> Vec<[f32; NVERTS]> {
        if let Some(p) = &self.params {
            cell_vcs
                .iter()
                .map(|vcs| {
                    let mut x_chemoas = [0.0f32; NVERTS];
                    vcs.iter().zip(x_chemoas.iter_mut()).for_each(
                        |(vc, x)| {
                            let r = p.center_mag
                                * p.slope
                                * (vc - &p.center).mag();
                            if r > 0.0 {
                                *x = r
                            }
                        },
                    );
                    x_chemoas
                })
                .collect()
        } else {
            vec![[0.0f32; NVERTS]; cell_vcs.len()]
        }
    }
}

#[derive(Clone)]
pub struct BdryGenerator {
    params: Option<BdryParams>,
}

impl BdryGenerator {
    pub fn new(params: Option<BdryParams>) -> BdryGenerator {
        BdryGenerator { params }
    }

    pub fn generate(
        &self,
        cell_vcs: &[[V2d; NVERTS]],
    ) -> Vec<[f32; NVERTS]> {
        if let Some(p) = &self.params {
            cell_vcs
                .iter()
                .map(|vcs| {
                    let mut x_bdrys = [0.0f32; NVERTS];
                    vcs.iter().zip(x_bdrys.iter_mut()).for_each(
                        |(vc, x)| {
                            if (p.skip_bb_check
                                && is_point_in_poly_no_bb_check(
                                    vc, &p.shape,
                                ))
                                || is_point_in_poly(
                                    vc, &p.bbox, &p.shape,
                                )
                            {
                                *x = p.mag;
                            }
                        },
                    );
                    x_bdrys
                })
                .collect()
        } else {
            vec![[0.0f32; NVERTS]; cell_vcs.len()]
        }
    }
}
