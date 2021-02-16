# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 22:10:15 2015

@author: brian
"""

import copy

import numpy as np

import environment
import geometry

output_mech_labels = [
    "x",
    "y",
    "edge_lengths",
    "sum_forces_x",
    "sum_forces_y",
    "edge_forces_plus_x",
    "edge_forces_plus_y",
    "edge_forces_minus_x",
    "edge_forces_minus_y",
    "rgtp_forces_x",
    "rgtp_forces_y",
    "cyto_forces_x",
    "cyto_forces_y",
    "local_strains",
    "uiv_x",
    "uiv_y",
]

output_chem_labels = [
    "rac_act",
    "rac_inact",
    "rac_cyto",
    "rho_act",
    "rho_inact",
    "rho_cyto",
    "old_x_coa",
    "old_x_cil",
    "kgtp_rac",
    "kgtp_rho",
    "kdgtp_rac",
    "kdgtp_rho",
    "rac_rands",
]

output_info_labels = output_mech_labels + output_chem_labels

for index, label in enumerate(output_info_labels):
    exec("{}_ix = {}".format(label, repr(index)))

num_info_labels = len(output_mech_labels + output_chem_labels)

info_indices_dict = {
    x: i for i, x in enumerate(output_mech_labels + output_chem_labels)
}

parameter_labels = ["cell_r",
                    "rest_edge_len",
                    "rest_area",
                    "stiffness_edge",
                    "const_protrusive",
                    "const_retractive",
                    "stiffness_cyto",
                    "k_mem_on_vertex",
                    "k_mem_off",
                    "diffusion_rgtp",
                    "init_rac",
                    "init_rho",
                    "halfmax_vertex_rgtp_act",
                    "halfmax_vertex_rgtp_conc",
                    "tot_rac",
                    "tot_rho",
                    "kgtp_rac",
                    "kgtp_rac_auto",
                    "kdgtp_rac",
                    "kdgtp_rho_on_rac",
                    "halfmax_tension_inhib",
                    "tension_inhib",
                    "kgtp_rho",
                    "kgtp_rho_auto",
                    "kdgtp_rho",
                    "kdgtp_rac_on_rho",
                    "randomization",
                    "rand_avg_t",
                    "rand_std_t",
                    "rand_mag",
                    "num_rand_vs",
                    "total_rgtp"]


# -----------------------------------------------------------------
def calc_coa_distrib_exp(coa_halfmax_dist):
    if coa_halfmax_dist > 1e-3:
        return np.log(0.5) / coa_halfmax_dist
    else:
        return 0.0


def refine_raw_params(raw_params):
    params = copy.deepcopy(raw_params)

    params["coa_mag"] = raw_params["coa_mag"] / 16
    params["coa_halfmax_dist"] = raw_params["coa_halfmax_dist"] / params["l"]
    params["coa_distrib_exp"] = calc_coa_distrib_exp(params["coa_halfmax_dist"])
    params["init_cyto_rgtp"] = 1 - raw_params["init_inact_rgtp"] - \
                               raw_params["init_act_rgtp"]
    # --------------
    params["kgtp_rac"] = (
            params["kgtp"] * raw_params["kgtp_rac"] * params["t"]
    )  # per second
    params["kdgtp_rac"] = (
            params["kgtp"] * raw_params["kdgtp_rac"] * params["t"]
    )  # per second
    # --------------
    params["kgtp_rho"] = (params["kgtp"] * raw_params["kgtp_rho"]) * params["t"]
    params["kdgtp_rho"] = (params["kgtp"] * raw_params["kdgtp_rho"]) * params[
        "t"]
    # --------------
    params["kgtp_rac_auto"] = \
        raw_params["kgtp_rac_autoact"] * params["kgtp"] * params["t"]
    params["kgtp_rho_auto"] = \
        raw_params["kgtp_rho_autoact"] * params["kgtp"] * params["t"]
    # --------------
    params["kdgtp_rho_on_rac"] = (
            params["kgtp"]
            * raw_params["kdgtp_rho_on_rac"] * params["t"]
    )  # per second
    params["kdgtp_rac_on_rho"] = (
            params["kgtp"]
            * raw_params["kdgtp_rac_on_rho"] * params["t"]
    )  # per second
    # --------------
    params["k_mem_off"] = (
            params["k_mem_off"] * raw_params["kgdi"] * params["t"]
    )  # per second
    params["k_mem_on_vertex"] = (
            params["k_mem_on"] * raw_params["kdgdi"] * params["t"] / 16
    )  # per second
    del params["k_mem_on"]
    # --------------
    params["diffusion_rgtp"] = \
        raw_params["diffusion_rgtp"] * (params["t"] / (params["l"] ** 2))
    # --------------
    params["halfmax_vertex_rgtp_act"] = \
        raw_params["halfmax_vertex_rgtp_act"] / 16
    # --------------
    params["cell_r"] = raw_params["cell_r"] / params["l"]
    vert_thetas = np.pi * np.linspace(0, 2, endpoint=False, num=16)
    cell_verts = np.transpose(
        np.array(
            [
                params["cell_r"] * np.cos(vert_thetas),
                params["cell_r"] * np.sin(vert_thetas),
            ]
        )
    )
    edge_vectors = geometry.calculate_edge_vectors(cell_verts)
    edge_lengths = geometry.calculate_vec_mags(edge_vectors)

    params["rest_edge_len"] = np.average(edge_lengths)
    params["rest_area"] = geometry.calculate_polygon_area(cell_verts)
    params["halfmax_vertex_rgtp_conc"] = params["halfmax_vertex_rgtp_act"] / \
                                         params["rest_edge_len"]

    params["vertex_eta"] = (params["vertex_eta"] * params["eta"] / 16) / \
                           (params["f"] / (params["l"] / params["t"]))
    params["stiffness_edge"] = \
        raw_params["stiffness_edge"] * \
        params["l3d"] * \
        params["l"] / params["f"]
    params["const_protrusive"] = \
        raw_params["lm_ss"] * raw_params["lm_h"] * \
        (params["rest_edge_len"] * params["l"]) * \
        params["halfmax_rgtp_max_f_frac"] / params["f"]
    params["const_retractive"] = (
            raw_params["force_rho"]
            * params["const_protrusive"]
    )
    params["stiffness_cyto"] = raw_params["stiffness_cyto"] / params["f"] / 16
    # --------------
    params["close_zero_at"] = raw_params["close_zero_at"] / params["l"]
    params["close_one_at"] = raw_params["close_one_at"] / params["l"]
    # --------------
    params["rand_avg_t"] = raw_params[
                               "rand_avg_t"
                           ] * 60.0 / params["t"]
    params["rand_std_t"] = raw_params[
                               "rand_std_t"
                           ] * 60.0 / params["t"]
    params["rand_mag"] = raw_params[
        "rand_mag"
    ]

    params["snap_period"] = raw_params["snap_period"] / params["t"]
    return params


# ==============================================================


def expand_x_cils_array(
        num_cell_groups, cell_group_defs, this_cell_group_def
):
    x_cils_def = this_cell_group_def["x_cils"]

    num_defs = len(
        list(
            x_cils_def.keys()))

    if num_defs != num_cell_groups:
        raise Exception(
            "Number of cell groups does not equal number of keys in x_cils_def."
        )

    x_cils = []
    for cgi in range(num_cell_groups):
        cg = cell_group_defs[cgi]
        cg_name = cg["group_name"]
        intercellular_contact_factor_mag = x_cils_def[
            cg_name]

        x_cils += (
                      cell_group_defs[cgi]["num_cells"]
                  ) * [intercellular_contact_factor_mag]

    return np.array(x_cils)
