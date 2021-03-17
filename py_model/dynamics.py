# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:22:43 2015

@author: Brian
"""

import hardio as fw
import utilities as general_utilities
import numba as nb
import numpy as np
import copy
import chemistry
import geometry
import mechanics
from hardio import BASIC_INFO, DYNAMICS_INFO, VOL_EX_INFO, DETAILED_VOL_EX_INFO


# ----------------------------------------------------------------------------------------
def pack_state_array(
        phase_var_indices,
        ode_cellwide_phase_var_indices,
        system_info_at_tstep):
    phase_var_array = (
        np.transpose(system_info_at_tstep[:, phase_var_indices])
    ).flatten()
    ode_cellwide_phase_var_array = system_info_at_tstep[
        0, ode_cellwide_phase_var_indices
    ]

    return np.append(phase_var_array, ode_cellwide_phase_var_array)


# -----------------------------------------------------------------------------------------
def interpret_state_array(rac_act_ix, rac_inact_ix, rho_act_ix,
                          rho_inact_ix, x_ix, y_ix, state_array):
    phase_vars = state_array

    rac_mem_active_start_ix = rac_act_ix * 16
    rac_mem_active_end_ix = rac_mem_active_start_ix + 16

    rac_acts = phase_vars[
               rac_mem_active_start_ix:rac_mem_active_end_ix
               ]

    rac_mem_inactive_start_ix = rac_inact_ix * 16
    rac_mem_inactive_end_ix = rac_mem_inactive_start_ix + 16

    rac_inacts = phase_vars[
                 rac_mem_inactive_start_ix:rac_mem_inactive_end_ix
                 ]

    rho_mem_active_start_ix = rho_act_ix * 16
    rho_mem_active_end_ix = rho_mem_active_start_ix + 16

    rho_acts = phase_vars[
               rho_mem_active_start_ix:rho_mem_active_end_ix
               ]

    rho_mem_inactive_start_ix = rho_inact_ix * 16
    rho_mem_inactive_end_ix = rho_mem_inactive_start_ix + 16

    rho_inacts = phase_vars[
                 rho_mem_inactive_start_ix:rho_mem_inactive_end_ix
                 ]

    x_start_ix = x_ix * 16
    x_end_ix = x_start_ix + 16

    x = phase_vars[x_start_ix:x_end_ix]

    y_start_ix = y_ix * 16
    y_end_ix = y_start_ix + 16

    y = phase_vars[y_start_ix:y_end_ix]

    poly = general_utilities.make_verts_array_given_xs_and_ys(x, y
                                                              )
    return rac_acts, rho_acts, rac_inacts, rho_inacts, poly


# ----------------------------------------------------------------------------------------
def unpack_state_array(num_phase_var_indices, state_array):
    # reversing append
    node_phase_var_array = state_array
    ode_cellwide_phase_vars = np.array([])

    # reversing flatten
    phase_vars = np.transpose(np.array(np.split(
        node_phase_var_array,
        num_phase_var_indices)))

    return phase_vars, ode_cellwide_phase_vars


# ----------------------------------------------------------------------------------------
def pack_state_array_from_system_history(
        phase_var_indices,
        ode_cellwide_phase_var_indices,
        system_info):
    state_array = pack_state_array(
        phase_var_indices,
        ode_cellwide_phase_var_indices,
        system_info)

    return state_array


# ----------------------------------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_sum(num_elements, sequence):
    result = 0
    for i in range(num_elements):
        result += sequence[i]

    return result


# ----------------------------------------------------------------------------------------
def eulerint(f, current_state, t0, t1, args, num_int_steps, cell_ix,
             curr_tpoint, rac_act_ix, rac_inact_ix,
             rho_act_ix, rho_inact_ix,
             x_ix, y_ix, writer):
    focus_verts = [0, 15]
    states = np.zeros(
        (2,
         current_state.shape[0]),
        dtype=np.float64)

    states[0] = copy.deepcopy(current_state)
    # logging.log(level=BASIC_INFO, msg="-----------------------------------")
    # logging.log(level=BASIC_INFO, msg="curr_tpoint: {}, cell: {}"
    #             .format(curr_tpoint, cell_ix))
    dt = (t1 - t0) / num_int_steps

    for int_step in range(num_int_steps):
        # logging.log(level=DYNAMICS_INFO,
        #             msg="-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+")
        _, _, _, _, init_poly \
            = interpret_state_array(rac_act_ix, rac_inact_ix,
                                    rho_act_ix, rho_inact_ix,
                                    x_ix, y_ix, current_state)
        # for ix in focus_verts:
        #     logging.log(level=DYNAMICS_INFO,
        #                 msg="init poly[{}]: {}".format(ix,
        #                                                init_poly[
        #                                                    ix]))
        deltas, sum_forces, edge_forces_plus, edge_forces_minus, \
        rgtp_forces, cyto_forces, verts_before_ve, verts_after_ve = f(
            focus_verts,
            curr_tpoint, int_step, dt, writer, cell_ix, current_state, *args)
        current_state = current_state + dt * deltas
        _, _, _, _, final_poly = interpret_state_array(rac_act_ix,
                                                       rac_inact_ix,
                                                       rho_act_ix,
                                                       rho_inact_ix,
                                                       x_ix, y_ix,
                                                       current_state)

        # for ix in focus_verts:
        #     logging \
        #         .log(level=DYNAMICS_INFO,
        #              msg="actual Delta poly({}): {}"
        #              .format(ix,
        #                      verts_before_ve[ix] -
        #                      init_poly[ix])
        #              )
        #     logging \
        #         .log(level=DYNAMICS_INFO,
        #              msg="Delta poly after VE({}): {}"
        #              .format(ix, verts_after_ve[ix] - init_poly[ix])
        #              )
        #     logging \
        #         .log(level=DYNAMICS_INFO,
        #              msg="final poly[{}]: {}"
        #              .format(ix, final_poly[0])
        #              )
        curr_tpoint += dt

    states[1] = copy.deepcopy(current_state)

    return states


def cell_dynamics(
        focus_verts,
        tpoint,
        int_step,
        dt,
        writer,
        cell_ix,
        state_array,
        num_cells,
        all_cells_verts,
        num_phase_vars,
        rac_act_ix,
        rest_edge_len,
        rac_inact_ix,
        rho_act_ix,
        rho_inact_ix,
        x_ix,
        y_ix,
        kgtp_rac,
        kdgtp_rac,
        kgtp_rho,
        kdgtp_rho,
        kgtp_rac_auto,
        kgtp_rho_auto,
        kdgtp_rho_on_rac,
        kdgtp_rac_on_rho,
        k_mem_on_vertex,
        k_mem_off,
        halfmax_vertex_rgtp_conc,
        diffusion_rgtp,
        vertex_eta,
        stiffness_edge,
        halfmax_vertex_rgtp_act,
        const_protrusive,
        const_retractive,
        rest_area,
        stiffness_cyto,
        x_coas,
        close_point_smoothness_factors,
        x_cils,
        halfmax_tension_inhib,
        tension_inhib,
        rac_rands,
        coa_updates,
        cil_updates
):
    phase_vars = state_array

    rac_mem_active_start_ix = rac_act_ix * 16
    rac_mem_active_end_ix = rac_mem_active_start_ix + 16

    rac_acts = phase_vars[
               rac_mem_active_start_ix:rac_mem_active_end_ix
               ]

    rac_mem_inactive_start_ix = rac_inact_ix * 16
    rac_mem_inactive_end_ix = rac_mem_inactive_start_ix + 16

    rac_inacts = phase_vars[
                 rac_mem_inactive_start_ix:rac_mem_inactive_end_ix
                 ]

    rho_mem_active_start_ix = rho_act_ix * 16
    rho_mem_active_end_ix = rho_mem_active_start_ix + 16

    rho_acts = phase_vars[
               rho_mem_active_start_ix:rho_mem_active_end_ix
               ]

    rho_mem_inactive_start_ix = rho_inact_ix * 16
    rho_mem_inactive_end_ix = rho_mem_inactive_start_ix + 16

    rho_inacts = phase_vars[
                 rho_mem_inactive_start_ix:rho_mem_inactive_end_ix
                 ]

    x_start_ix = x_ix * 16
    x_end_ix = x_start_ix + 16

    x = phase_vars[x_start_ix:x_end_ix]

    y_start_ix = y_ix * 16
    y_end_ix = y_start_ix + 16

    y = phase_vars[y_start_ix:y_end_ix]

    poly = general_utilities.make_verts_array_given_xs_and_ys(x, y
                                                              )

    rac_cyto = (
            1
            - calculate_sum(16, rac_acts)
            - calculate_sum(16, rac_inacts)
    )
    rho_cyto = (
            1
            - calculate_sum(16, rho_acts)
            - calculate_sum(16, rho_inacts)
    )

    sum_forces, edge_forces_plus, edge_forces_minus, uevs, rgtp_forces, \
    cyto_forces, \
    edge_strains, local_strains, \
    uivs = mechanics.calculate_forces(
        poly,
        rac_acts,
        rho_acts,
        rest_edge_len,
        stiffness_edge,
        halfmax_vertex_rgtp_act,
        const_protrusive,
        const_retractive,
        rest_area,
        stiffness_cyto,
    )

    sum_forces_x = sum_forces[:, 0]
    sum_forces_y = sum_forces[:, 1]

    only_tensile_local_strains = np.zeros_like(local_strains)
    for i in range(16):
        local_strain = local_strains[i]
        if local_strain > 0:
            only_tensile_local_strains[i] = local_strain

    edgeplus_lengths = geometry.calculate_edgeplus_lengths(poly)
    avg_edge_lengths = geometry.calculate_average_edge_length_around_nodes(
        edgeplus_lengths
    )
    conc_rac_acts = chemistry.calc_concs(rac_acts, avg_edge_lengths)

    kgtps_rac = chemistry.calculate_kgtp_rac(
        conc_rac_acts,
        halfmax_vertex_rgtp_conc,
        kgtp_rac,
        kgtp_rac_auto,
        x_coas,
        rac_rands,
        x_cils,
        close_point_smoothness_factors,
    )

    conc_rho_acts = chemistry.calc_concs(rho_acts, avg_edge_lengths
                                         )

    global_tension = np.sum(only_tensile_local_strains) / 16
    if global_tension < 0.0:
        global_tension = 0.0
    strain_inhibition = tension_inhib * \
                        chemistry.hill_function3(
                            halfmax_tension_inhib,
                            global_tension
                        )

    kdgtps_rac = chemistry.calculate_kdgtp_rac(
        conc_rho_acts,
        halfmax_vertex_rgtp_conc,
        kdgtp_rac,
        kdgtp_rho_on_rac,
        x_cils,
        halfmax_tension_inhib,
        tension_inhib,
        only_tensile_local_strains,
    )

    kgtps_rho = chemistry.calculate_kgtp_rho(
        conc_rho_acts,
        x_cils,
        halfmax_vertex_rgtp_conc,
        kgtp_rho,
        kgtp_rho_auto,
    )

    kdgtps_rho = chemistry.calculate_kdgtp_rho(
        conc_rac_acts,
        halfmax_vertex_rgtp_conc,
        kdgtp_rho,
        kdgtp_rac_on_rho,
    )

    conc_rac_inacts = chemistry.calc_concs(rac_inacts, avg_edge_lengths)
    conc_rho_inact = chemistry.calc_concs(rho_inacts, avg_edge_lengths)

    rac_act_net_fluxes = chemistry.calculate_net_fluxes(
        conc_rac_acts,
        diffusion_rgtp,
        edgeplus_lengths,
    )
    rac_inact_net_fluxes = chemistry.calculate_net_fluxes(
        conc_rac_inacts,
        diffusion_rgtp,
        edgeplus_lengths,
    )
    rho_act_net_fluxes = chemistry.calculate_net_fluxes(
        conc_rho_acts,
        diffusion_rgtp,
        edgeplus_lengths,
    )
    rho_inact_net_fluxes = chemistry.calculate_net_fluxes(
        conc_rho_inact,
        diffusion_rgtp,
        edgeplus_lengths,
    )

    delta_rac_activated = np.zeros(16, dtype=np.float64)
    delta_rac_inactivated = np.zeros(16, dtype=np.float64)

    delta_rac_cytosol_to_membrane = np.zeros(16, dtype=np.float64)

    delta_rho_activated = np.zeros(16, dtype=np.float64)
    delta_rho_inactivated = np.zeros(16, dtype=np.float64)

    delta_rho_cytosol_to_membrane = np.zeros(16, dtype=np.float64)

    delta_x = np.zeros(16, dtype=np.float64)
    delta_y = np.zeros(16, dtype=np.float64)
    new_verts = np.zeros((16, 2), dtype=np.float64)
    np.zeros(2, dtype=np.float64)
    np.zeros(2, dtype=np.float64)

    # logging.log(level=DYNAMICS_INFO,
    #             msg="tstep: {}, int_step: {}".format(tpoint, int_step))
    # logging.log(level=DYNAMICS_INFO, msg="eta: {}".format(vertex_eta))
    # logging.log(level=DYNAMICS_INFO, msg="1/eta: {}".format(1 / vertex_eta))
    # for ix in focus_verts:
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="rgtp_forces[{}]: {}".format(ix, rgtp_forces[ix]))
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="edge_forces[{}]: {}".format(ix,
    #                                                  edge_forces_plus[ix]))
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="cyto_forces[{}]: {}".format(ix, cyto_forces[ix]))
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="expected sum forces ({}) = {}".format(ix,
    #                                                            rgtp_forces[
    #                                                                ix] +
    #                                                            edge_forces_plus[
    #                                                                ix] +
    #                                                            edge_forces_minus[
    #                                                                ix] +
    #                                                            cyto_forces[
    #                                                                ix]))
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="sum_forces[{}]: {}".format(ix, sum_forces[ix]))

    poly_area = geometry.calculate_polygon_area(poly)
    data = [("tpoint", tpoint),
            ("poly", [[float(v) for v in x] for x in poly]),
            ("rac_acts", [float(v) for v in rac_acts]),
            ("rac_inacts", [float(v) for v in rac_inacts]),
            ("rho_acts", [float(v) for v in rho_acts]),
            ("rho_inacts", [float(v) for v in rho_inacts]),
            ("sum_forces", [list([float(x), float(y)]) for x, y in
                            zip(sum_forces_x, sum_forces_y)]),
            ("uivs", [[float(v) for v in x] for x in uivs]),
            ("rgtp_forces", [[float(v) for v in x] for x in rgtp_forces]),
            ("edge_forces", [[float(v) for v in x] for x in edge_forces_plus]),
            ("cyto_forces", [[float(v) for v in x] for x in cyto_forces]),
            ("kgtps_rac", [float(v) for v in kgtps_rac]),
            ("kdgtps_rac", [float(v) for v in kdgtps_rac]),
            ("kgtps_rho", [float(v) for v in kgtps_rho]),
            ("kdgtps_rho", [float(v) for v in kdgtps_rho]),
            ("x_cils", [float(v) for v in x_cils]),
            ("x_coas", [float(v) for v in x_coas]),
            ("edge_strains", [float(v) for v in local_strains]),
            ("avg_tens_strain", [float(global_tension) for _ in local_strains]),
            ("poly_area", poly_area),
            ("rac_act_net_fluxes", [float(v) for v in rac_act_net_fluxes]),
            ("rac_inact_net_fluxes", [float(v) for v in rac_inact_net_fluxes]),
            ("rho_act_net_fluxes", [float(v) for v in rho_act_net_fluxes]),
            ("rho_inact_net_fluxes", [float(v) for v in rho_inact_net_fluxes]),
            ("x_tens", [float(strain_inhibition) for _ in local_strains])]
    # for d in data:
    #     logging.log(level=99, msg="{}: {}".format(d[0], d[1]))
    writer.save_int_step(data)

    for ni in range(16):
        old_coord = poly[ni]

        new_verts[ni][0] = old_coord[0] + dt * sum_forces_x[ni] / vertex_eta
        new_verts[ni][1] = old_coord[1] + dt * sum_forces_y[ni] / vertex_eta

    # for ix in focus_verts:
    #     logging.log(level=DYNAMICS_INFO, msg="delta.poly[{}]: {}"
    #                 .format(ix, (new_verts[0] - poly[0]) / dt)
    #                 )
    #     logging.log(level=DYNAMICS_INFO,
    #                 msg="expected Delta poly 0 ({}): {}"
    #                 .format(ix, dt * sum_forces[0] / vertex_eta)
    #                 )

    verts_before_ve = copy.deepcopy(new_verts)

    # calculate volume exclusion effects
    num_bisection_iterations = 10
    max_movement_mag = dt * const_protrusive / vertex_eta

    for other_ci in range(num_cells):
        if other_ci != cell_ix:
            # logging.log(level=VOL_EX_INFO, msg="testing poly: {}".format(
            #     other_ci))
            # logging.log(level=VOL_EX_INFO,
            #             msg="coords: {}".format(all_cells_verts[other_ci]))
            # logging.log(level=VOL_EX_INFO,
            #             msg="max movement mag: {}".format(max_movement_mag))
            are_new_nodes_inside_other_cell = \
                geometry.are_points_inside_polygon(
                    new_verts, all_cells_verts[other_ci]
                )
            # logging.log(level=VOL_EX_INFO, msg="in poly: {}".format(
            #     [i for (i, x) in
            #      enumerate(are_new_nodes_inside_other_cell) if x])
            #             )
            for ni in range(16):
                if are_new_nodes_inside_other_cell[ni]:
                    # logging.log(level=VOL_EX_INFO,
                    #             msg="fixing vertex {} violation (current: {})".format(
                    #                 ni, new_verts[ni]))
                    new_verts[ni] = enforce_volume_exclusion_for_vertex(
                        poly[ni],
                        new_verts[ni],
                        uivs[ni],
                        all_cells_verts[other_ci],
                        num_bisection_iterations,
                        max_movement_mag,
                    )
    verts_after_ve = copy.deepcopy(new_verts)

    for ni in range(16):
        new_coord = new_verts[ni]
        old_coord = poly[ni]

        delta_x[ni] = (new_coord[0] - old_coord[0]) / dt
        delta_y[ni] = (new_coord[1] - old_coord[1]) / dt

    for ni in range(16):
        # finish assigning chemistry variables
        delta_rac_activated[ni] = kgtps_rac[ni] * rac_inacts[ni]
        delta_rac_inactivated[ni] = kdgtps_rac[ni] * rac_acts[ni]

        delta_rac_on = k_mem_on_vertex * rac_cyto
        delta_rac_off = k_mem_off * rac_inacts[ni]
        delta_rac_cytosol_to_membrane[ni] = delta_rac_on - delta_rac_off

        delta_rho_activated[ni] = kgtps_rho[ni] * rho_inacts[ni]
        delta_rho_inactivated[ni] = kdgtps_rho[ni] * rho_acts[ni]

        delta_rho_on = k_mem_on_vertex * rho_cyto
        delta_rho_off = k_mem_off * rho_inacts[ni]
        delta_rho_cytosol_to_membrane[ni] = delta_rho_on - delta_rho_off

    # set up ode array
    ode_array = np.empty(num_phase_vars * 16)

    for i in range(16):
        ode_array[i] = (
                delta_rac_activated[i]
                - delta_rac_inactivated[i]
                + rac_act_net_fluxes[i]
        )

        ode_array[i + 16] = (
                delta_rac_inactivated[i]
                - delta_rac_activated[i]
                + rac_inact_net_fluxes[i]
                + delta_rac_cytosol_to_membrane[i]
        )

        ode_array[i + 2 * 16] = (
                delta_rho_activated[i]
                - delta_rho_inactivated[i]
                + rho_act_net_fluxes[i]
        )

        ode_array[i + 3 * 16] = (
                delta_rho_inactivated[i]
                - delta_rho_activated[i]
                + rho_inact_net_fluxes[i]
                + delta_rho_cytosol_to_membrane[i]
        )

        ode_array[i + 4 * 16] = delta_x[i]

        ode_array[i + 5 * 16] = delta_y[i]

    return ode_array, sum_forces, edge_forces_plus, edge_forces_minus, \
           rgtp_forces, cyto_forces, verts_before_ve, verts_after_ve


# -----------------------------------------------------------------
def enforce_volume_exclusion_for_vertex(
        old_coord,
        new_coord,
        unit_inside_pointing_vector,
        polygon,
        num_bisection_iterations,
        max_movement_mag,
):
    # min_x, max_x, min_y, max_y = geometry.calculate_polygon_bb(
    #     polygon)

    is_old_in_poly = geometry.is_point_in_polygon_without_bb_check(
        old_coord, polygon
    )

    while is_old_in_poly:
        old_coord = old_coord + max_movement_mag * \
                    unit_inside_pointing_vector

        # logging.log(level=DETAILED_VOL_EX_INFO, msg="trial old v: {}".format(
        #     old_coord))
        # num_bisection_iterations = int(num_bisection_iterations*1.5)
        is_old_in_poly = geometry.is_point_in_polygon_without_bb_check(
            old_coord, polygon)

    # if we have reached here, then we know that the old_coord is in the
    # polygon, and the new coord is not in the polygon
    ok_coord = old_coord

    # logging.log(level=DETAILED_VOL_EX_INFO,
    #             msg="settling with okay v: {} (in poly: {})".format(
    #                 old_coord,
    #                 geometry.is_point_in_polygon_without_bb_check(
    #                     old_coord,
    #                     polygon)))
    problem_coord = new_coord
    np.zeros(2, dtype=np.float64)

    # logging.log(level=DETAILED_VOL_EX_INFO,
    #             msg="problem v: {}".format(problem_coord))
    for i in range(num_bisection_iterations):
        test_coord = 0.5 * (ok_coord + problem_coord)

        # logging.log(level=DETAILED_VOL_EX_INFO,
        #             msg="testing: {}".format(test_coord))

        if geometry.is_point_in_polygon_without_bb_check(
                test_coord, polygon
        ):
            # logging.log(level=DETAILED_VOL_EX_INFO, msg="setting as problem")
            problem_coord = test_coord
        else:
            # logging.log(level=DETAILED_VOL_EX_INFO, msg="setting as ok")
            ok_coord = test_coord

    # logging.log(level=DETAILED_VOL_EX_INFO,
    #             msg="returning ok: {}".format(ok_coord))
    return ok_coord
