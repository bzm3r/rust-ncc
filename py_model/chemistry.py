# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:27:54 2015
@author: Brian
"""

import numba as nb
import numpy as np
from hardio import COA_INFO


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def hill_function3(thresh, sig):
    pow_sig = sig ** 3.0
    pow_thresh = thresh ** 3.0

    return pow_sig / (pow_thresh + pow_sig)


# -----------------------------------------------------------------
# @nb.jit(nopython=True)


def calculate_kgtp_rac(
        rac_act_net_fluxes,
        halfmax_vertex_rgtp_conc,
        kgtp_rac,
        kgtp_rac_auto,
        x_coas,
        randomization_factors,
        intercellular_contact_factors,
        close_point_smoothness_factors,
):
    num_vertices = rac_act_net_fluxes.shape[0]
    result = np.empty(num_vertices, dtype=np.float64)

    for i in range(num_vertices):
        i_plus1 = (i + 1) % num_vertices
        i_minus1 = (i - 1) % num_vertices

        cil_factor = (
                             intercellular_contact_factors[i]
                             + intercellular_contact_factors[i_plus1]
                             + intercellular_contact_factors[i_minus1]
                     ) / 3.0
        smooth_factor = np.max(close_point_smoothness_factors[i])
        x_coa = x_coas[i]

        if cil_factor > 0.0 or smooth_factor > 1e-6:
            x_coa = 0.0

        rac_autoact_hill_effect = hill_function3(
            halfmax_vertex_rgtp_conc,
            rac_act_net_fluxes[i])
        kgtp_rac_autoact = (
                kgtp_rac_auto
                * rac_autoact_hill_effect
        )

        if kgtp_rac_autoact > 1.25 * kgtp_rac_auto:
            kgtp_rac_autoact = 1.25 * kgtp_rac_auto

        kgtp_rac_base = (
                                1.0 + randomization_factors[
                            i] + x_coa) * kgtp_rac
        result[i] = kgtp_rac_base + kgtp_rac_autoact

    return result


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_kgtp_rho(
        conc_rho_act,
        intercellular_contact_factors,
        halfmax_vertex_rgtp_conc,
        kgtp_rho,
        kgtp_rho_auto,
):
    result = np.empty(16)
    for i in range(16):
        kgtp_rho_autoact = kgtp_rho_auto * hill_function3(
            halfmax_vertex_rgtp_conc,
            conc_rho_act[i]
        )

        i_plus1 = (i + 1) % 16
        i_minus1 = (i - 1) % 16

        cil_factor = (
                             intercellular_contact_factors[i]
                             + intercellular_contact_factors[i_plus1]
                             + intercellular_contact_factors[i_minus1]
                     ) / 3.0

        result[i] = (
                            1.0 + cil_factor
                    ) * kgtp_rho + kgtp_rho_autoact
    return result


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_kdgtp_rac(
        conc_rho_acts,
        halfmax_vertex_rgtp_conc,
        kdgtp_rac,
        kdgtp_rho_on_rac,
        intercellular_contact_factors,
        halfmax_tension_inhib,
        tension_inhib,
        local_strains,
):
    result = np.empty(16, dtype=np.float64)

    global_tension = np.sum(local_strains) / 16
    if global_tension < 0.0:
        global_tension = 0.0
    strain_inhibition = tension_inhib * \
                        hill_function3(
                            halfmax_tension_inhib,
                            global_tension
                        )

    for i in range(16):
        kdgtp_rho_mediated_rac_inhib = (
                kdgtp_rho_on_rac
                * hill_function3(
            halfmax_vertex_rgtp_conc,
            conc_rho_acts[i],
        )
        )

        i_plus1 = (i + 1) % 16
        i_minus1 = (i - 1) % 16

        cil_factor = (
                             intercellular_contact_factors[i]
                             + intercellular_contact_factors[i_plus1]
                             + intercellular_contact_factors[i_minus1]
                     ) / 3.0

        result[i] = (
                            1.0 + cil_factor + strain_inhibition
                    ) * kdgtp_rac + kdgtp_rho_mediated_rac_inhib

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_kdgtp_rho(
        rac_act_net_fluxes,
        halfmax_vertex_rgtp_conc,
        kdgtp_rho,
        kdgtp_rac_on_rho,
):
    result = np.empty(16, dtype=np.float64)

    for i in range(16):
        kdgtp_rac_mediated_rho_inhib = (
                kdgtp_rac_on_rho
                * hill_function3(
            halfmax_vertex_rgtp_conc,
            rac_act_net_fluxes[i],
        )
        )

        result[i] = kdgtp_rho + kdgtp_rac_mediated_rho_inhib

    return result


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calc_concs(species, avg_edge_lengths):
    result = np.empty(16, dtype=np.float64)

    for i in range(16):
        result[i] = species[i] / avg_edge_lengths[i]

    return result


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_flux_terms(concentrations, diffusion_rgtp, edgeplus_lengths):
    result = np.empty(16, dtype=np.float64)

    for i in range(16):
        i_plus1_ix = (i + 1) % 16

        result[i] = (
                - diffusion_rgtp
                * (concentrations[i_plus1_ix] - concentrations[i])
                / edgeplus_lengths[i]
        )

    return result


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_diffusion(concentrations, diffusion_rgtp, edgeplus_lengths):
    result = np.empty(16, dtype=np.float64)

    fluxes = calculate_flux_terms(
        concentrations,
        diffusion_rgtp,
        edgeplus_lengths,
    )

    for i in range(16):
        i_minus1_ix = (i - 1) % 16

        result[i] = fluxes[i_minus1_ix] - fluxes[i]

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calc_x_cils(
        this_cell_ix,
        num_cells,
        cil_mag,
        close_point_smoothness_factors,
):
    x_cils = np.zeros(16, dtype=np.float64)

    for other_ci in range(num_cells):
        if other_ci != this_cell_ix:
            for ni in range(16):
                x_cils[ni] = (
                        cil_mag
                        * close_point_smoothness_factors[ni][other_ci]
                )

    return x_cils


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calculate_x_coas(
        num_cells,
        this_cell_ix,
        coa_distrib_exp,
        coa_mag,
        intercellular_dist_squared_matrix,
        line_segment_intersection_matrix,
        coa_los_penalty,
):
    x_coas = np.zeros(16, dtype=np.float64)
    too_close_dist_squared = 1e-6

    for ni in range(16):
        this_node_x_coa = x_coas[ni]

        this_node_relevant_line_seg_intersection_slice = \
            line_segment_intersection_matrix[
                ni]
        this_node_relevant_dist_squared_slice = \
            intercellular_dist_squared_matrix[ni]

        for other_ci in range(num_cells):
            if other_ci != this_cell_ix:
                this_node_other_cell_relevant_line_seg_intersection_slice = \
                    this_node_relevant_line_seg_intersection_slice[
                        other_ci]
                relevant_dist_squared_slice = \
                    this_node_relevant_dist_squared_slice[
                        other_ci]
                for other_ni in range(16):
                    line_segment_between_node_intersects_polygon = \
                        this_node_other_cell_relevant_line_seg_intersection_slice[
                            other_ni]
                    intersection_factor = (
                            1.0
                            / (
                                    line_segment_between_node_intersects_polygon + 1.0)
                            ** coa_los_penalty
                    )

                    dist_squared_between_nodes = \
                        relevant_dist_squared_slice[other_ni]

                    # logging.log(level=COA_INFO, msg="====================")
                    # logging.log(level=COA_INFO,
                    #             msg="(ci: {}, vi: {}, ovi: {}, oci: {}):".format(
                    #                 this_cell_ix, ni, other_ci, other_ni))
                    # logging.log(level=COA_INFO, msg="dist: {}".format(np.sqrt(
                    #     dist_squared_between_nodes)))
                    # logging.log(level=COA_INFO, msg="num_intersects: {}".format(
                    #     line_segment_between_node_intersects_polygon))
                    # logging.log(level=COA_INFO,
                    #             msg="coa_los_penalty: {}".format(
                    #                 intersection_factor))

                    coa_signal = 0.0
                    if dist_squared_between_nodes > too_close_dist_squared:
                        coa_signal = (
                                np.exp(
                                    coa_distrib_exp
                                    * np.sqrt(dist_squared_between_nodes)
                                )
                                * intersection_factor
                        )
                    old_x_coa = this_node_x_coa
                    this_node_x_coa += coa_mag * coa_signal
                    # logging.log(level=COA_INFO,
                    #             msg="coa_signal: {}".format(coa_signal))
                    # logging.log(level=COA_INFO,
                    #             msg="coa_mag: {}".format(coa_mag))
                    # logging.log(level=COA_INFO,
                    #             msg="new = {} + {}".format(old_x_coa,
                    #                                        coa_mag * coa_signal))

        x_coas[ni] = this_node_x_coa

    return x_coas


# -------------------------------------------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_chemoattractant_shielding_effect_factors(
        this_cell_ix,
        num_cells,
        intercellular_dist_squared_matrix,
        line_segment_intersection_matrix,
        chemoattractant_shielding_effect_length,
):
    chemoattractant_shielding_effect_factors = np.zeros(
        16, dtype=np.float64)

    for ni in range(16):
        this_node_relevant_line_seg_intersection_slice = \
            line_segment_intersection_matrix[
                ni]
        this_node_relevant_dist_squared_slice = \
            intercellular_dist_squared_matrix[ni]

        sum_weights = 0.0

        for other_ci in range(num_cells):
            if other_ci != this_cell_ix:
                this_node_other_cell_relevant_line_seg_intersection_slice = \
                    this_node_relevant_line_seg_intersection_slice[
                        other_ci]
                this_node_other_cell_relevant_dist_squared_slice = \
                    this_node_relevant_dist_squared_slice[
                        other_ci]

                for other_ni in range(16):
                    if (
                            this_node_other_cell_relevant_line_seg_intersection_slice[
                                other_ni] == 0):
                        ds = this_node_other_cell_relevant_dist_squared_slice[
                            other_ni]
                        sum_weights += np.exp(
                            np.log(0.25)
                            * (ds / chemoattractant_shielding_effect_length)
                        )

        chemoattractant_shielding_effect_factors[ni] = 1.0 / \
                                                       (1.0 + sum_weights)

    return chemoattractant_shielding_effect_factors
