// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::{P2D, calc_poly_area, max_f32};
use crate::utils::{circ_ix_plus, circ_ix_minus};

pub fn calc_edge_unit_vecs(vertex_coords: &[P2D]) -> Vec<P2D> {
    let nvs = vertex_coords.len();
    (0..nvs)
        .map(|i| {
            let plus_i = circ_ix_plus(i, nvs);
            (vertex_coords[plus_i] - vertex_coords[i]).unitize()
        })
        .collect()
}

pub fn calc_edge_lens(edge_unit_vecs: &[P2D]) -> Vec<f32> {
    edge_unit_vecs.iter().map(|euv| euv.mag()).collect()
}

pub fn calc_edge_strains(edge_lens: &[f32], rel: f32) -> Vec<f32> {
    edge_lens.iter().map(|&el| (el / rel) - 1.0).collect()
}

pub fn calc_global_strain(edge_lens: &[f32], rel: f32, nverts: usize) -> f32 {
    edge_lens.iter().sum::<f32>() / (nverts as f32 * rel)
}

pub fn calc_edge_forces(
    edge_strains: &[f32],
    edge_unit_vecs: &[P2D],
    stiffness_edge: f32,
) -> Vec<P2D> {
    let nvs = edge_strains.len();
    (0..nvs)
        .map(|i| {
            edge_unit_vecs[i].scalar_mul(edge_strains[i] * stiffness_edge)
        })
        .collect()
}

pub fn calc_cyto_forces(vertex_coords: &[P2D], unit_inward_vecs: &[P2D], rest_area: f32, stiffness_cyto: f32) -> Vec<P2D> {
    let area = calc_poly_area(vertex_coords);
    let areal_strain = area/rest_area - 1.0;
    let mag = stiffness_cyto * areal_strain / (vertex_coords.len() as f32);
    unit_inward_vecs.iter().map(|uiv| uiv.scalar_mul(mag)).collect()
}

pub fn calc_rgtp_forces(
    rac_acts: &[f32],
    rho_acts: &[f32],
    unit_inward_vecs: &[P2D],
    max_protrusive_f: f32,
    max_retractive_f: f32,
    max_f_activity: f32,
) -> Vec<P2D> {
    let const_protrusive = max_protrusive_f / max_f_activity;
    let const_retractive = max_retractive_f / max_f_activity;
    unit_inward_vecs
        .iter()
        .zip(rac_acts.iter().zip(rho_acts.iter()))
        .map(|(uiv, (ra, pa))| {
            if ra > pa {
                uiv * -1.0 * max_f32(ra - pa, max_f_activity) * const_protrusive
            } else {
                uiv * max_f32(pa - ra, max_f_activity) * const_retractive
            }
        })
        .collect()
}


