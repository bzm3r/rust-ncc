// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

mod gen_coa;
// pub mod dat2d;
pub mod dat_4d;
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
    PhysContactFactors, PhysicalContactGenerator,
};
use crate::interactions::RelativeRgtpActivity::{
    RacDominant, RhoDominant,
};
use crate::math::geometry::{BBox, Poly};
use crate::math::v2d::V2D;
use crate::parameters::InteractionParams;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

/// The relative Rho GTPase activity at a cell is positive if Rac1
/// dominates, otherwise it is negative.
#[derive(Copy, Clone, Debug, Deserialize, Serialize)]
pub enum RelativeRgtpActivity {
    RhoDominant(f64),
    RacDominant(f64),
}

impl RelativeRgtpActivity {
    pub fn to_f64(&self) -> f64 {
        match self {
            RelativeRgtpActivity::RacDominant(value) => *value,
            RelativeRgtpActivity::RhoDominant(value) => *value,
        }
    }

    pub fn from_f64(value: f64) -> Self {
        if value > 0.0 {
            RacDominant(value)
        } else {
            RhoDominant(-1.0 * value)
        }
    }

    /// Edges are always defined as going from the smaller indexed
    /// vertex to the bigger indexed vertex.
    pub fn mix_rel_rgtp_act_across_edge(
        smaller_vertex: RelativeRgtpActivity,
        bigger_vertex: RelativeRgtpActivity,
        t: f64,
    ) -> RelativeRgtpActivity {
        RelativeRgtpActivity::from_f64(
            smaller_vertex.to_f64() * (1.0 - t)
                + bigger_vertex.to_f64() * t,
        )
    }
}

#[derive(
    Copy, Clone, Debug, Default, Deserialize, Serialize, PartialEq,
)]
pub struct Interactions {
    pub x_cals: [f64; NVERTS],
    pub x_cils: [f64; NVERTS],
    pub x_adhs: [V2D; NVERTS],
    pub x_chem_attrs: [f64; NVERTS],
    pub x_coas: [f64; NVERTS],
    pub x_bdrys: [f64; NVERTS],
}

/// Generates interaction related factors.
#[derive(Clone, Deserialize, Serialize)]
pub struct InteractionGenerator {
    /// Vertex coordinates, per cell, for all cells in the simulation.
    cell_polys: Vec<Poly>,
    all_rgtps: Vec<[RelativeRgtpActivity; NVERTS]>,
    /// Generates CIL/CAL related interaction information. In other
    /// words, interactions that require cells to engage in physical
    /// contact.
    pub phys_contact_generator: PhysicalContactGenerator,
    coa_generator: Option<CoaGenerator>,
    chem_attr_generator: Option<ChemAttrGenerator>,
    bdry_generator: Option<BdryEffectGenerator>,
}

pub struct ContactData {
    pub oci: usize,
    pub poly: Poly,
}

impl InteractionGenerator {
    pub fn new(
        cell_verts: &[[V2D; NVERTS]],
        cell_rgtps: &[[RelativeRgtpActivity; NVERTS]],
        params: InteractionParams,
    ) -> InteractionGenerator {
        let cell_polys = cell_verts
            .iter()
            .map(|vs| Poly::from_points(vs))
            .collect::<Vec<Poly>>();
        let phys_contact_generator = PhysicalContactGenerator::new(
            &cell_polys,
            params.phys_contact,
        );
        let coa_generator = params.coa.map(|coa_params| {
            CoaGenerator::new(&cell_polys, coa_params)
        });
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
        self.cell_polys[cell_ix] = Poly::from_points(vs);
        if let Some(coa_gen) = self.coa_generator.as_mut() {
            coa_gen.update(cell_ix, &self.cell_polys)
        }
        if let Some(_chema_gen) = self.chem_attr_generator.as_mut() {
            unimplemented!()
        }
        if let Some(_bdry_gen) = self.bdry_generator.as_mut() {
            unimplemented!()
        }
        self.phys_contact_generator
            .update(cell_ix, &self.cell_polys);
    }

    pub fn generate(
        &self,
        rel_rgtps: &Vec<[RelativeRgtpActivity; NVERTS]>,
    ) -> Vec<Interactions> {
        let num_cells = self.cell_polys.len();
        let PhysContactFactors { adh, cil, cal } =
            self.phys_contact_generator.generate(rel_rgtps);
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
            .map(|ci| Interactions {
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
                poly: self.cell_polys[oci],
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
