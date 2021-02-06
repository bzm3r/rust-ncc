// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::capped_linear_fn;
use crate::math::geometry::calc_poly_area;
use crate::math::v2d::V2D;
use crate::utils::circ_ix_plus;
use crate::NVERTS;

/// Calculate edge vectors of a polygon.
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

/// Calculate elastic forces due to stretching of edges.
pub fn calc_edge_forces(
    edge_strains: &[f64; NVERTS],
    edge_unit_vecs: &[V2D; NVERTS],
    stiffness_edge: f64,
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    (0..NVERTS).for_each(|i| {
        // Elastic relationship: stiffness * strain = magnitude of
        // force. Direction of force is along the unite edge vector.
        r[i] = edge_strains[i] * stiffness_edge * edge_unit_vecs[i]
    });
    r
}

/// Calculate forces generated due to deformation of the cytoplasm.
/// Assume cytoplasm is a "elastic" liquid, where the force generated
/// is linearly proportional to the areal strain:
/// `(final_area - initial_area)/initial_area`
pub fn calc_cyto_forces(
    vertex_coords: &[V2D; NVERTS],
    unit_inward_vecs: &[V2D; NVERTS],
    rest_area: f64,
    stiffness_cyto: f64,
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    let area = calc_poly_area(vertex_coords);
    let areal_strain = (area / rest_area) - 1.0;
    let mag = stiffness_cyto * areal_strain;
    (0..NVERTS).for_each(|i| r[i] = mag * unit_inward_vecs[i]);
    r
}

/// Except for `X_acts`, and `unit_inward_vecs`, the rest of the
/// inputs to this function come from parameters. We need the
/// active fraction of a Rho GTPase at a vertex to determine the force
/// exerted due to signalling of this fraction on downstream
/// force-generating elements. For simplicity, we assume a piece-wise
/// linear (approximating a sigmoid) relationship between Rho GTPase
/// activity, and the force generated. The shape of the sigmoid is
/// governed by `halfmax_vertex_rgtp_act`.
pub fn calc_rgtp_forces(
    rac_acts: &[f64; NVERTS],
    rho_acts: &[f64; NVERTS],
    unit_inward_vecs: &[V2D; NVERTS],
    halfmax_vertex_rgtp_act: f64,
    const_protrusive: f64,
    const_retractive: f64,
) -> [V2D; NVERTS] {
    let mut r = [V2D::default(); NVERTS];
    for i in 0..NVERTS {
        // Direction force will point in.
        let uiv = unit_inward_vecs[i];
        let ra = rac_acts[i];
        let pa = rho_acts[i];

        // If Rac1 activity is greater than RhoA activity, then
        // generate a protrusive force; one that points in the
        // direction `-1.0 * uiv`. Else, generate a contractile force
        // which points in the direction `uiv`.
        let mag = if ra > pa {
            -1.0 * const_protrusive
                * capped_linear_fn(
                    ra - pa,
                    0.0,
                    2.0 * halfmax_vertex_rgtp_act,
                )
        } else {
            const_retractive
                * capped_linear_fn(
                    pa - ra,
                    0.0,
                    2.0 * halfmax_vertex_rgtp_act,
                )
        };

        r[i] = mag * uiv;
    }
    r
}
