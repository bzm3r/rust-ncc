// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

mod gen_coa;
// pub mod dat2d;
pub mod dat4d;
pub mod dat_sym2d;
pub mod dat_sym4d;
mod dat_utils;
pub mod gen_bdry;
pub mod gen_chemoa;
mod gen_phys;

use crate::interactions::dat_sym2d::SymCcDat;
use crate::interactions::gen_bdry::BdryEffectGenerator;
use crate::interactions::gen_chemoa::ChemAttrGenerator;
use crate::interactions::gen_coa::CoaGenerator;
use crate::interactions::gen_phys::{
    ClosePoint, PhysContactFactors, PhysicalContactGenerator,
};
use crate::math::geometry::BBox;
use crate::math::v2d::V2D;
use crate::parameters::InteractionParams;
use crate::NVERTS;
use avro_schema_derive::Schematize;
use serde::{Deserialize, Serialize};

pub type RgtpState = f32;

#[derive(
    Copy, Clone, Debug, Default, Deserialize, Schematize, Serialize,
)]
pub struct CellInteractions {
    pub x_cals: [f32; NVERTS],
    pub x_cils: [f32; NVERTS],
    pub x_adhs: [V2D; NVERTS],
    pub x_chem_attrs: [f32; NVERTS],
    pub x_coas: [f32; NVERTS],
    pub x_bdrys: [f32; NVERTS],
}

/// Generates interaction related factors.
#[derive(Clone)]
pub struct InteractionGenerator {
    /// Vertex coordinates, per cell, for all cells in the simulation.
    cell_polys: Vec<[V2D; NVERTS]>,
    all_rgtps: Vec<[RgtpState; NVERTS]>,
    /// Generates CIL/CAL related interaction information. In other words,
    /// interactions that require cells to engage in physical contact.
    phys_contact_generator: PhysicalContactGenerator,
    coa_generator: Option<CoaGenerator>,
    chem_attr_generator: Option<ChemAttrGenerator>,
    bdry_generator: Option<BdryEffectGenerator>,
}

pub struct ContactData {
    pub oci: usize,
    pub close_verts: Vec<usize>,
    pub poly: [V2D; NVERTS],
    pub poly_bb: BBox,
}

impl InteractionGenerator {
    pub fn new(
        cell_polys: &[[V2D; NVERTS]],
        cell_rgtps: &[[RgtpState; NVERTS]],
        params: InteractionParams,
    ) -> InteractionGenerator {
        let cell_bbs = cell_polys
            .iter()
            .map(|vs| BBox::from_points(vs))
            .collect::<Vec<BBox>>();
        let phys_contact_generator = PhysicalContactGenerator::new(
            &cell_bbs,
            cell_polys,
            params.phys_contact,
        );
        let coa_generator = params
            .coa
            .map(|p| CoaGenerator::new(&cell_bbs, cell_polys, p));
        let chem_attr_generator =
            params.chem_attr.map(ChemAttrGenerator::new);
        let bdry_generator =
            params.bdry.map(BdryEffectGenerator::new);
        InteractionGenerator {
            cell_polys: cell_polys.iter().copied().collect(),
            all_rgtps: cell_rgtps.iter().copied().collect(),
            phys_contact_generator,
            coa_generator,
            chem_attr_generator,
            bdry_generator,
        }
    }

    pub fn update(&mut self, cell_ix: usize, vs: &[V2D; NVERTS]) {
        self.cell_polys[cell_ix] = *vs;
        let bb = BBox::from_points(vs);
        if let Some(coa_gen) = self.coa_generator.as_mut() {
            coa_gen.update(cell_ix, &bb, vs, &self.cell_polys)
        }
        if let Some(_chema_gen) = self.chem_attr_generator.as_mut() {
            unimplemented!()
        }
        if let Some(_bdry_gen) = self.bdry_generator.as_mut() {
            unimplemented!()
        }
        self.phys_contact_generator.update(
            cell_ix,
            &bb,
            vs,
            &self.cell_polys,
        );
    }

    pub fn generate(&self) -> Vec<CellInteractions> {
        let num_cells = self.cell_polys.len();
        let PhysContactFactors { adh, cil, cal } =
            self.phys_contact_generator.generate(&self.all_rgtps);
        let r_coas = self
            .coa_generator
            .as_ref()
            .map_or(vec![[0.0; NVERTS]; num_cells], |gen| {
                gen.generate()
            });
        let r_chemoas = self
            .chem_attr_generator
            .as_ref()
            .map_or(vec![[0.0; NVERTS]; num_cells], |gen| {
                gen.generate(&self.cell_polys)
            });
        let r_bdrys = self
            .bdry_generator
            .as_ref()
            .map_or(vec![[0.0; NVERTS]; num_cells], |gen| {
                gen.generate(&self.cell_polys)
            });
        (0..num_cells)
            .map(|ci| CellInteractions {
                x_cals: cal[ci],
                x_cils: cil[ci],
                x_adhs: adh[ci],
                x_chem_attrs: r_chemoas[ci],
                x_coas: r_coas[ci],
                x_bdrys: r_bdrys[ci],
            })
            .collect()
    }

    pub fn get_physical_contacts(&self, ci: usize) -> Vec<usize> {
        let num_cells = self.cell_polys.len();
        (0..num_cells)
            .filter(|&oci| {
                self.phys_contact_generator.contacts.get(ci, oci)
            })
            .collect()
    }

    pub fn get_contact_data(&self, ci: usize) -> Vec<ContactData> {
        self.get_physical_contacts(ci)
            .into_iter()
            .map(|oci| ContactData {
                oci,
                close_verts: self
                    .phys_contact_generator
                    .get_close_verts(ci, oci),
                poly: self.cell_polys[oci],
                poly_bb: self.phys_contact_generator.contact_bbs[oci],
            })
            .collect()
    }
}

/// Generate a `SymCcDat<bool>` which records whether two cells
/// are roughly in contact, based on intersection of their contact
/// bounding boxes.
pub fn generate_contacts(contact_bbs: &[BBox]) -> SymCcDat<bool> {
    let mut contacts = SymCcDat::new(contact_bbs.len(), false);
    for (ci, bb) in contact_bbs.iter().enumerate() {
        for (oxi, obb) in contact_bbs[(ci + 1)..].iter().enumerate() {
            contacts.set(ci, ci + 1 + oxi, obb.intersects(bb));
        }
    }
    contacts
}
