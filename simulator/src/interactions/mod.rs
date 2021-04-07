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
use crate::interactions::gen_phys::PhysicalContactGenerator;
use crate::interactions::RelativeRgtpActivity::{
    RacDominant, RhoDominant,
};
use crate::math::geometry::{BBox, Poly};
use crate::math::v2d::V2d;
use crate::parameters::InteractionParams;
use crate::NVERTS;
use serde::{Deserialize, Serialize};

/// The relative Rho GTPase activity at a cell is positive if Rac1
/// dominates, otherwise it is negative.
#[derive(Copy, Clone, Debug, Deserialize, Serialize, PartialEq)]
pub enum RelativeRgtpActivity {
    RhoDominant(f64),
    RacDominant(f64),
}

impl Default for RelativeRgtpActivity {
    fn default() -> Self {
        RelativeRgtpActivity::RhoDominant(0.0)
    }
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
            RhoDominant(value)
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
    pub x_adhs: [V2d; NVERTS],
    pub x_chem_attrs: [f64; NVERTS],
    pub x_coas: [f64; NVERTS],
    pub x_bdrys: [f64; NVERTS],
}

/// Generates interaction related factors.
#[derive(Clone, Deserialize, Serialize)]
pub struct InteractionGenerator {
    all_rgtps: Vec<[RelativeRgtpActivity; NVERTS]>,
    /// Generates CIL/CAL related interaction information. In other
    /// words, interactions that require cells to engage in physical
    /// contact.
    pub phys_contact_generator: PhysicalContactGenerator,
    coa_generator: Option<CoaGenerator>,
    chem_attr_generator: Option<ChemAttrGenerator>,
    bdry_generator: Option<BdryEffectGenerator>,
}

pub struct Contact {
    pub oci: usize,
    pub poly: Poly,
}

impl InteractionGenerator {
    pub fn new(
        cell_verts: &[[V2d; NVERTS]],
        cell_rgtps: &[[RelativeRgtpActivity; NVERTS]],
        params: InteractionParams,
    ) -> InteractionGenerator {
        let cell_polys = cell_verts
            .iter()
            .map(|vs| Poly::from_verts(vs))
            .collect::<Vec<Poly>>();
        let phys_contact_generator = PhysicalContactGenerator::new(
            &cell_polys,
            cell_rgtps,
            params.phys_contact,
        );
        let coa_generator = params.coa.map(|coa_params| {
            CoaGenerator::new(
                &cell_polys,
                coa_params,
                &phys_contact_generator,
            )
        });
        let chem_attr_generator =
            params.chem_attr.map(ChemAttrGenerator::new);
        let bdry_generator =
            params.bdry.map(BdryEffectGenerator::new);
        InteractionGenerator {
            all_rgtps: cell_rgtps.iter().copied().collect(),
            phys_contact_generator,
            coa_generator,
            chem_attr_generator,
            bdry_generator,
        }
    }

    pub fn update_phys_for(
        &mut self,
        cell_ix: usize,
        cell_poly: Poly,
        rel_rgtps: [RelativeRgtpActivity; NVERTS],
    ) {
        self.phys_contact_generator
            .update_for(cell_ix, cell_poly, rel_rgtps);
    }

    pub fn update_coa_for(
        &mut self,
        cell_ix: usize,
        cell_poly: Poly,
    ) {
        if let Some(coa_gen) = self.coa_generator.as_mut() {
            coa_gen.update_for(
                cell_ix,
                cell_poly,
                &self.phys_contact_generator,
            )
        }
    }

    pub fn update_chemoa_for(
        &mut self,
        _cell_ix: usize,
        _cell_poly: Poly,
    ) {
        unimplemented!()
    }

    pub fn get_physical_contacts(&self, ci: usize) -> Vec<usize> {
        let num_cells = self.phys_contact_generator.cell_polys.len();
        (0..num_cells)
            .filter(|&oci| {
                self.phys_contact_generator
                    .contact_matrix
                    .get(ci, oci)
            })
            .collect()
    }

    pub fn get_contacts(&self, ci: usize) -> Vec<Contact> {
        self.get_physical_contacts(ci)
            .into_iter()
            .map(|oci| Contact {
                oci,
                poly: self.phys_contact_generator.cell_polys[oci],
            })
            .collect()
    }
}

/// Generate a `SymCcDat<bool>` which records whether two cells
/// are roughly in contact, based on intersection of their contact
/// bounding boxes.
pub fn gen_contact_matrix(contact_bbs: &[BBox]) -> SymCcDat<bool> {
    let mut contacts = SymCcDat::new(contact_bbs.len(), false);
    for (ci, bb) in contact_bbs.iter().enumerate() {
        for (oxi, obb) in contact_bbs[(ci + 1)..].iter().enumerate() {
            contacts.set(ci, ci + 1 + oxi, obb.intersects(bb));
        }
    }
    contacts
}
