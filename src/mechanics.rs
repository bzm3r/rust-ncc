// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::consts::NVERTS;
use crate::math::{calc_poly_area, capped_linear_function, P2D};
use crate::utils::circ_ix_plus;

pub fn calc_edge_vecs(vertex_coords: &[P2D; NVERTS as usize]) -> [P2D; NVERTS as usize] {
    let mut r = [P2D::default(); NVERTS as usize];
    (0..NVERTS as usize).for_each(|i| {
        let plus_i = circ_ix_plus(i, NVERTS as usize);
        r[i] = vertex_coords[plus_i] - vertex_coords[i];
    });
    r
}

pub fn calc_edge_forces(
    edge_strains: &[f32; NVERTS as usize],
    edge_unit_vecs: &[P2D; NVERTS as usize],
    stiffness_edge: f32,
) -> [P2D; NVERTS as usize] {
    let mut r = [P2D::default(); NVERTS as usize];
    (0..NVERTS as usize)
        .for_each(|i| r[i] = edge_strains[i] * stiffness_edge * edge_unit_vecs[i]);
    r
}

pub fn calc_cyto_forces(
    vertex_coords: &[P2D; NVERTS as usize],
    unit_inward_vecs: &[P2D; NVERTS as usize],
    rest_area: f32,
    stiffness_cyto: f32,
) -> [P2D; NVERTS as usize] {
    let area = calc_poly_area(vertex_coords);
    let areal_strain = (area / rest_area) - 1.0;
    let mag = stiffness_cyto * areal_strain / (vertex_coords.len() as f32);
    let mut r = [P2D::default(); NVERTS as usize];
    (0..NVERTS as usize).for_each(|i| r[i] = mag * unit_inward_vecs[i]);
    r
}

pub fn calc_rgtp_forces(
    rac_acts: &[f32; NVERTS as usize],
    rho_acts: &[f32; NVERTS as usize],
    unit_inward_vecs: &[P2D; NVERTS as usize],
    halfmax_vertex_rgtp_act: f32,
    const_protrusive: f32,
    const_retractive: f32,
) -> [P2D; NVERTS as usize] {
    let nvs = unit_inward_vecs.len();
    let mut r = [P2D::default(); NVERTS as usize];
    for i in 0..nvs {
        let uiv = unit_inward_vecs[i];
        let ra = rac_acts[i];
        let pa = rho_acts[i];

        let mag = if ra > pa {
            -1.0 * const_protrusive * capped_linear_function(ra - pa, 2.0 * halfmax_vertex_rgtp_act)
        } else {
            const_retractive * capped_linear_function(pa - ra, 2.0 * halfmax_vertex_rgtp_act)
        };

        r[i] = mag * uiv;
    }
    r
}
