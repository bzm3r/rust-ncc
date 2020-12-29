// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::capped_linear_function;
use crate::math::geometry::calc_poly_area;
use crate::math::v2d::V2D;
use crate::utils::circ_ix_plus;
use crate::NVERTS;

pub fn calc_edge_vecs(
    vertex_coords: &[V2D; NVERTS],
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    (0..NVERTS).for_each(|i| {
        let plus_i = circ_ix_plus(i, NVERTS);
        r[i] = vertex_coords[plus_i] - vertex_coords[i];
    });
    r
}

pub fn calc_edge_forces(
    edge_strains: &[f32; NVERTS],
    edge_unit_vecs: &[V2D; NVERTS],
    stiffness_edge: f32,
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    (0..NVERTS).for_each(|i| {
        r[i] = edge_strains[i] * stiffness_edge * edge_unit_vecs[i]
    });
    r
}

pub fn calc_cyto_forces(
    vertex_coords: &[V2D; NVERTS],
    unit_inward_vecs: &[V2D; NVERTS],
    rest_area: f32,
    stiffness_cyto: f32,
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    let area = calc_poly_area(vertex_coords);
    let areal_strain = (area / rest_area) - 1.0;
    let mag = stiffness_cyto * areal_strain;
    (0..NVERTS).for_each(|i| r[i] = mag * unit_inward_vecs[i]);
    r
}

pub fn calc_rgtp_forces(
    rac_acts: &[f32; NVERTS],
    rho_acts: &[f32; NVERTS],
    unit_inward_vecs: &[V2D; NVERTS],
    halfmax_vertex_rgtp_act: f32,
    const_protrusive: f32,
    const_retractive: f32,
) -> [V2D; NVERTS] {
    let nvs = unit_inward_vecs.len();
    let mut r = [V2D::default(); NVERTS];
    for i in 0..nvs {
        let uiv = unit_inward_vecs[i];
        let ra = rac_acts[i];
        let pa = rho_acts[i];

        let mag = if ra > pa {
            -1.0 * const_protrusive
                * capped_linear_function(
                    ra - pa,
                    2.0 * halfmax_vertex_rgtp_act,
                )
        } else {
            const_retractive
                * capped_linear_function(
                    pa - ra,
                    2.0 * halfmax_vertex_rgtp_act,
                )
        };

        r[i] = mag * uiv;
    }
    r
}
