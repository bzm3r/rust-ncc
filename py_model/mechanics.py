# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:26:37 2015

@author: Brian
"""

import numba as nb
import numpy as np

import geometry


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def capped_linear_function(max_x, x):
    if x > max_x:
        return 1.0
    else:
        return x / max_x


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def linear_function(max_x, x):
    return x / max_x


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_phys_space_bdry_contact_factors(
        this_cell_coords, space_physical_bdry_polygons
):
    if space_physical_bdry_polygons.size == 0:
        return np.zeros(len(this_cell_coords))
    else:
        return geometry.are_points_inside_polygons(
            this_cell_coords, space_physical_bdry_polygons
        )


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_cytoplasmic_force(
        this_cell_coords,
        rest_area,
        stiffness_cyto,
        uivs,
):
    current_area = abs(geometry.calculate_polygon_area(this_cell_coords))

    area_strain = (current_area - rest_area) / rest_area
    force_mag = area_strain * stiffness_cyto

    return geometry.multiply_vectors_by_scalar(
        uivs, force_mag
    )


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_spring_edge_forces(
        this_cell_coords, stiffness_edge, rest_edge_len
):
    edge_vectors_to_plus = np.empty((16, 2), dtype=np.float64)
    #edge_vectors_to_minus = np.empty((16, 2), dtype=np.float64)

    for i in range(16):
        i_plus_1 = (i + 1) % 16
        i_minus_1 = (i - 1) % 16

        edge_vectors_to_plus[i] = \
            this_cell_coords[i_plus_1] - this_cell_coords[i]
        # edge_vectors_to_minus[i] = \
        #     this_cell_coords[i_minus_1] - this_cell_coords[i]

    plus_dirn_edge_length = np.linalg.norm(edge_vectors_to_plus, axis=1)

    # minus_dirn_edge_length = np.linalg.norm(edge_vectors_to_minus, axis=1)

    edge_strains_plus = np.empty(16, dtype=np.float64)
    # edge_strains_minus = np.empty(16, dtype=np.float64)
    local_average_strains = np.empty(16, dtype=np.float64)

    for i in range(16):
        edge_strain_plus = (
                                   plus_dirn_edge_length[i] - rest_edge_len
                           ) / rest_edge_len
        # edge_strain_minus = (
        #                             minus_dirn_edge_length[i] - rest_edge_len
        #                     ) / rest_edge_len

        edge_strains_plus[i] = edge_strain_plus
        # edge_strains_minus[i] = edge_strain_minus

    edge_strains_minus = np.roll(edge_strains_plus, shift=1, axis=0)
    local_average_strains = (edge_strains_plus + edge_strains_minus) * 0.5
    unit_edge_disp_vecs_plus = edge_vectors_to_plus / plus_dirn_edge_length[
                                                      :, np.newaxis]
    # unit_edge_disp_vecs_minus = edge_vectors_to_minus / minus_dirn_edge_length[
    #                                                   :, np.newaxis]

    edge_forces_plus_mags = np.zeros(16, dtype=np.float64)
    # edge_forces_minus_mags = np.zeros(16, dtype=np.float64)
    for i in range(16):
        edge_forces_plus_mags[i] = edge_strains_plus[i] * stiffness_edge
        #edge_forces_minus_mags[i] = edge_strains_minus[i] * stiffness_edge

    edge_forces_plus = edge_forces_plus_mags[:, np.newaxis] * unit_edge_disp_vecs_plus
    edge_forces_minus = -1.0 * np.roll(edge_forces_plus, shift=1, axis=0)

    #edge_forces_minus = edge_forces_minus_mags[:, np.newaxis] *
    # unit_edge_disp_vecs_minus

    #norm_sum_edge_fs = np.linalg.norm(np.sum(
    # edge_forces_plus, axis=0) +
    #                                  np.sum(
    #                                  edge_forces_minus, axis=0))
    #print("norm_sum_edge_fs: {}".format(norm_sum_edge_fs))

    return local_average_strains, edge_forces_plus, edge_forces_minus, \
           edge_strains_plus, unit_edge_disp_vecs_plus, edge_vectors_to_plus


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def determine_rac_rho_domination(rac_acts, rho_acts):
    domination_array = np.empty(16, dtype=np.int64)

    for ni in range(16):
        if rac_acts[ni] < rho_acts[ni]:
            domination_array[ni] = 0
        else:
            domination_array[ni] = 1

    return domination_array


# -----------------------------------------------------------------


# @nb.jit(nopython=True)
def calculate_rgtpase_mediated_forces(
        rac_acts,
        rho_acts,
        halfmax_vertex_rgtp,
        const_protrusive,
        const_retractive,
        uivs,
):
    rgtpase_mediated_force_mags = np.zeros(16, dtype=np.float64)

    for ni in range(16):
        rac_activity = rac_acts[ni]
        rho_activity = rho_acts[ni]

        if rac_activity > rho_activity:
            mag_frac = capped_linear_function(
                2 * halfmax_vertex_rgtp, rac_activity - rho_activity
            )
            force_mag = const_protrusive * mag_frac
        else:
            force_mag = (
                    -1
                    * const_retractive
                    * capped_linear_function(
                2 * halfmax_vertex_rgtp, rho_activity - rac_activity
            )
            )

        rgtpase_mediated_force_mags[ni] = -1 * force_mag

    np.empty((16, 2), dtype=np.float64)
    result = geometry.multiply_vectors_by_scalars(
        uivs, rgtpase_mediated_force_mags
    )

    return result


# ----------------------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_forces(
        this_cell_coords,
        rac_acts,
        rho_acts,
        rest_edge_len,
        stiffness_edge,
        halfmax_vertex_rgtp,
        const_protrusive,
        const_retractive,
        rest_area,
        stiffness_cyto,
):
    uivs = geometry.calculate_unit_inside_pointing_vecs(
        this_cell_coords)

    rgtpase_mediated_forces = calculate_rgtpase_mediated_forces(
        rac_acts,
        rho_acts,
        halfmax_vertex_rgtp,
        const_protrusive,
        const_retractive,
        uivs,
    )

    cyto_forces = calculate_cytoplasmic_force(
        this_cell_coords,
        rest_area,
        stiffness_cyto,
        uivs,
    )

    local_average_strains, edge_forces_plus, edge_forces_minus, \
    edge_strains, uevs, edge_vecs_plus = \
        calculate_spring_edge_forces(
        this_cell_coords, stiffness_edge, rest_edge_len
    )

    sum_forces = rgtpase_mediated_forces + edge_forces_plus + \
                 edge_forces_minus + cyto_forces

    return (
        sum_forces,
        edge_forces_plus,
        edge_forces_minus,
        uevs,
        rgtpase_mediated_forces,
        cyto_forces,
        edge_strains,
        local_average_strains,
        uivs,
    )

# =============================================================================
