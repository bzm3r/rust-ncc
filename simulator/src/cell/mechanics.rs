// Copyright © 2020 Brian Merchant.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::math::capped_linear_fn;
use crate::math::geometry::calc_poly_area;
use crate::math::v2d::V2d;
use crate::utils::circ_ix_plus;
use crate::NVERTS;

/// Calculate edge vectors of a polygon.
pub fn calc_edge_vecs(
    vertex_coords: &[V2d; NVERTS],
) -> [V2d; NVERTS] {
    let mut r = [V2d::default(); NVERTS];
    (0..NVERTS).for_each(|i| {
        let plus_i = circ_ix_plus(i, NVERTS);
        r[i] = vertex_coords[plus_i] - vertex_coords[i];
    });
    r
}

/// Calculate elastic forces due to stretching of edges.
pub fn calc_edge_forces(
    edge_strains: &[f64; NVERTS],
    unit_edge_vecs: &[V2d; NVERTS],
    stiffness_edge: f64,
) -> [V2d; NVERTS] {
    let mut r = [V2d::default(); NVERTS];
    (0..NVERTS).for_each(|i| {
        // Elastic relationship: stiffness * strain = magnitude of
        // force. Direction of force is along the unite edge vector.
        r[i] = edge_strains[i] * stiffness_edge * unit_edge_vecs[i]
    });
    r
}

/// Calculate forces generated due to deformation of the cytoplasm.
/// Assume cytoplasm is a "elastic" liquid, where the force generated
/// is linearly proportional to the areal strain:
/// `(final_area - initial_area)/initial_area`
pub fn calc_cyto_forces(
    vertex_coords: &[V2d; NVERTS],
    unit_inward_vecs: &[V2d; NVERTS],
    rest_area: f64,
    stiffness_cyto: f64,
) -> [V2d; NVERTS] {
    let mut r = [V2d::default(); NVERTS];
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
/// governed by `halfmax_vertex_rgtp`.
pub fn calc_rgtp_forces(
    rac_acts: &[f64; NVERTS],
    rho_acts: &[f64; NVERTS],
    unit_inward_vecs: &[V2d; NVERTS],
    halfmax_vertex_rgtp: f64,
    const_protrusive: f64,
    const_retractive: f64,
) -> [V2d; NVERTS] {
    let mut r = [V2d::default(); NVERTS];
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
                    2.0 * halfmax_vertex_rgtp,
                )
        } else {
            const_retractive
                * capped_linear_fn(
                    pa - ra,
                    0.0,
                    2.0 * halfmax_vertex_rgtp,
                )
        };

        r[i] = mag * uiv;
    }
    r
}
