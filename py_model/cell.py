import utilities as utils
import numba as nb
import numpy as np

import chemistry
import dynamics
import geometry
import mechanics
import parameters

import copy

"""
Cell.

"""


# =============================================

class Cell:
    def __init__(
            self,
            cell_group_ix,
            cell_ix,
            params,
    ):
        self.cell_group_ix = cell_group_ix
        self.cell_ix = cell_ix

        self.num_tsteps = params["num_tsteps"]
        self.num_tpoints = self.num_tsteps + 1
        self.num_int_steps = params["num_int_steps"]

        self.curr_tstep = 0
        self.num_cells = params["num_cells"]

        self.tot_rac = params["tot_rac"]
        self.tot_rho = params["tot_rho"]

        self.curr_state = np.zeros(
            (
                16,
                len(parameters.output_info_labels),
            )
        )

        self.l = params["l"]
        self.t = params["t"]
        self.f = params["f"]

        self.vertex_eta = params["vertex_eta"]
        self.curr_verts = params["init_verts"]
        self.radius_resting = params["cell_r"]
        self.rest_edge_len = params["rest_edge_len"]
        edgeplus_lengths = geometry.calculate_edgeplus_lengths(
            self.curr_verts)
        self.init_average_edge_lengths = np.average(
            geometry.calculate_average_edge_length_around_nodes(
                edgeplus_lengths))
        self.rest_area = params["rest_area"]
        # ======================================================

        self.close_zero_at = params["close_zero_at"] ** 2
        self.close_one_at = params["close_one_at"] ** 2

        self.kgtp_rac = params["kgtp_rac"]
        self.kgtp_rac_auto = params["kgtp_rac_auto"]

        self.halfmax_vertex_rgtp_act = params["halfmax_vertex_rgtp_act"]
        self.halfmax_vertex_rgtp_conc = params["halfmax_vertex_rgtp_conc"]

        self.kdgtp_rac = params["kdgtp_rac"]
        self.kdgtp_rho_on_rac = params["kdgtp_rho_on_rac"]

        self.halfmax_tension_inhib = params[
            "halfmax_tension_inhib"
        ]
        self.tension_inhib = params[
            "tension_inhib"
        ]
        self.k_mem_off = params["k_mem_off"]
        self.k_mem_on_vertex = params["k_mem_on_vertex"]

        self.diffusion_rgtp = params["diffusion_rgtp"]
        # ======================================================

        self.kgtp_rho = params["kgtp_rho"]
        self.kgtp_rho_auto = params["kgtp_rho_auto"]

        self.kdgtp_rho = params["kdgtp_rho"]
        self.kdgtp_rac_on_rho = params["kdgtp_rac_on_rho"]

        # ==============================================================
        self.cil_mag = params["cil_mag"]
        self.coa_mag = params["coa_mag"]

        self.coa_mag = params["coa_mag"]
        self.coa_distrib_exp = params["coa_distrib_exp"]
        self.coa_los_penalty = params[
            "coa_los_penalty"]
        # =============================================================

        self.stiffness_edge = params["stiffness_edge"]
        self.stiffness_cyto = params["stiffness_cyto"]
        self.const_protrusive = params["const_protrusive"]
        self.const_retractive = params["const_retractive"]
        # =============================================================

        self.phase_var_indices = [
            parameters.rac_act_ix,
            parameters.rac_inact_ix,
            parameters.rho_act_ix,
            parameters.rho_inact_ix,
            parameters.x_ix,
            parameters.y_ix,
        ]
        self.num_phase_vars = len(self.phase_var_indices)
        self.total_num_phase_vars = self.num_phase_vars * \
                                    16

        self.initialize_phase_var_indices()

        # =============================================================
        self.rand_avg_t = params["rand_avg_t"]
        self.rand_std_t = params["rand_std_t"]
        self.rand_mag = params["rand_mag"]
        self.num_rand_vs = params["num_rand_vs"]
        self.randomization = params["randomization"]
        self.rac_rands = self.renew_rac_rands(0)
        self.next_randomization_event_tpoint = 1200

        # =============================================================

        self.all_cellwide_phase_var_indices = [
            parameters.rac_cyto_ix,
            parameters.rho_cyto_ix,
        ]
        self.ode_cellwide_phase_var_indices = []
        self.num_all_cellwide_phase_vars = len(
            self.all_cellwide_phase_var_indices)
        self.num_ode_cellwide_phase_vars = len(
            self.ode_cellwide_phase_var_indices)

        self.initialize_all_cellwide_phase_var_indices()
        self.initialize_ode_cellwide_phase_var_indices()

        # =============================================================

        self.pars_indices = [
            parameters.kgtp_rac_ix,
            parameters.kgtp_rho_ix,
            parameters.kdgtp_rac_ix,
            parameters.kdgtp_rho_ix,
            parameters.local_strains_ix,
            parameters.old_x_cil_ix,
        ]

        self.initialize_pars_indices()

        # =============================================================
        self.init_inact_rgtp = params[
            "init_inact_rgtp"
        ]
        self.init_act_rgtp = params[
            "init_act_rgtp"
        ]
        self.rgtp_distrib_def_for_randomization = [
            "unbiased random",
            0.0,
            0.0,
        ]
        self.init_rac = params["init_rac"]
        self.init_rho = params["init_rho"]

    # -----------------------------------------------------------------
    def insert_state_array_into_system_history(self, state_array):
        phase_vars, ode_cellwide_phase_vars = dynamics.unpack_state_array(
            self.num_phase_vars, state_array)
        self.curr_state[:, self.phase_var_indices] = phase_vars
        self.curr_state[
            0, self.ode_cellwide_phase_var_indices] = ode_cellwide_phase_vars

    # -----------------------------------------------------------------
    def initialize_cell(
            self,
            close_point_smoothness_factors, x_cils, x_coas,
    ):
        verts = self.curr_verts

        self.curr_state[:, [parameters.x_ix, parameters.y_ix]] = verts
        self.curr_state[:, parameters.edge_lengths_ix] = \
            self.rest_edge_len * np.ones(16)

        if self.cell_ix == 0:
            self.set_rgtp_distrib(
                self.init_rac, self.init_rho,
            )
        else:
            self.set_rgtp_distrib(
                self.init_rho, self.init_rac,
            )


        rac_acts = self.curr_state[:, parameters.rac_act_ix]
        rho_acts = self.curr_state[:, parameters.rho_act_ix]

        self.curr_state[:, parameters.old_x_coa_ix] = x_coas
        self.curr_state[:, parameters.old_x_cil_ix] = x_cils

        sum_forces, edge_forces_plus, edge_forces_minus, uevs, rgtpase_forces, \
        cyto_forces, \
        edge_strains, local_strains, uivs = \
            mechanics.calculate_forces(
                verts,
                rac_acts,
                rho_acts,
                self.rest_edge_len,
                self.stiffness_edge,
                self.halfmax_vertex_rgtp_act,
                self.const_protrusive,
                self.const_retractive,
                self.rest_area,
                self.stiffness_cyto,
            )

        self.curr_state[:, parameters.local_strains_ix] = local_strains

        edgeplus_lengths = geometry.calculate_edgeplus_lengths(verts)
        avg_edge_lengths = geometry.calculate_average_edge_length_around_nodes(
            edgeplus_lengths
        )

        conc_rac_acts = chemistry.calc_concs(rac_acts,
                                             avg_edge_lengths)

        conc_rho_acts = chemistry.calc_concs(rho_acts,
                                             avg_edge_lengths)

        self.curr_state[:, parameters.rac_rands_ix] = self.rac_rands

        kgtp_racs = chemistry.calculate_kgtp_rac(
            conc_rac_acts,
            self.halfmax_vertex_rgtp_conc,
            self.kgtp_rac,
            self.kgtp_rac_auto,
            x_coas,
            self.rac_rands,
            x_cils,
            close_point_smoothness_factors,
        )

        self.curr_state[:, parameters.kgtp_rac_ix] = kgtp_racs

        kgtp_rhos = chemistry.calculate_kgtp_rho(
            conc_rho_acts,
            x_cils,
            self.halfmax_vertex_rgtp_conc,
            self.kgtp_rho,
            self.kgtp_rho_auto,
        )
        self.curr_state[:, parameters.kgtp_rho_ix] = kgtp_rhos

        self.curr_state[:, parameters.kdgtp_rac_ix] = \
            chemistry.calculate_kdgtp_rac(
                conc_rho_acts,
                self.halfmax_vertex_rgtp_conc,
                self.kdgtp_rac,
                self.kdgtp_rho_on_rac,
                x_cils,
                self.halfmax_tension_inhib,
                self.tension_inhib,
                np.array([ls if ls > 0 else 0.0 for ls in local_strains]),
            )

        self.curr_state[:, parameters.kdgtp_rho_ix] = \
            chemistry.calculate_kdgtp_rho(
                conc_rac_acts,
                self.halfmax_vertex_rgtp_conc,
                self.kdgtp_rho,
                self.kdgtp_rac_on_rho,
            )

        # update mechanics parameters
        self.curr_state[:, [parameters.sum_forces_x_ix,
                            parameters.sum_forces_y_ix]] \
            = sum_forces
        self.curr_state[:, [parameters.edge_forces_plus_x_ix,
                            parameters.edge_forces_plus_y_ix]] = \
            edge_forces_plus
        self.curr_state[:, [parameters.edge_forces_minus_x_ix,
                            parameters.edge_forces_minus_y_ix]] = \
            edge_forces_minus
        self.curr_state[:, [parameters.rgtp_forces_x_ix,
                            parameters.rgtp_forces_y_ix]] \
            = rgtpase_forces
        self.curr_state[:, [parameters.cyto_forces_x_ix,
                            parameters.cyto_forces_y_ix]] = \
            cyto_forces

        self.curr_state[:, [parameters.uiv_x_ix, parameters.uiv_y_ix]] = uivs

        self.curr_verts = verts

    # -----------------------------------------------------------------
    def initialize_phase_var_indices(self):
        for index, sys_info_ix in enumerate(self.phase_var_indices):
            label = parameters.output_info_labels[sys_info_ix]
            setattr(self, "" + label + "_ix", index)

    # -----------------------------------------------------------------
    def initialize_pars_indices(self):
        for index, sys_info_ix in enumerate(self.pars_indices):
            label = parameters.output_info_labels[sys_info_ix]
            setattr(self, "" + label + "_ix", index)

    # -----------------------------------------------------------------
    def initialize_ode_cellwide_phase_var_indices(self):
        for index, sys_info_ix in enumerate(
                self.ode_cellwide_phase_var_indices):
            label = parameters.output_info_labels[sys_info_ix]
            setattr(self, "cellwide_" + label + "_ix", index)

    # -----------------------------------------------------------------
    def initialize_all_cellwide_phase_var_indices(self):
        for index, sys_info_ix in enumerate(
                self.all_cellwide_phase_var_indices):
            label = parameters.output_info_labels[sys_info_ix]
            setattr(self, "cellwide_" + label + "_ix", index)

    # -----------------------------------------------------------------
    def calculate_when_randomization_event_occurs(self):
        return self.curr_tstep + 1200

    # -----------------------------------------------------------------
    def set_rgtp_distrib(
            self,
            init_rac,
            init_rho,
    ):
        self.curr_state[:, parameters.rac_act_ix] = init_rac
        self.curr_state[:, parameters.rac_inact_ix] = init_rac
        self.curr_state[:, parameters.rho_act_ix] = init_rho
        self.curr_state[:, parameters.rho_inact_ix] = init_rho
        self.curr_state[:, parameters.rac_cyto_ix] = 1.0 - 2 * init_rac
        self.curr_state[:, parameters.rho_cyto_ix] = 1.0 - 2 * init_rho

    # -----------------------------------------------------------------

    def renew_rac_rands(self, tstep):
        if self.randomization:
            possible_rfs = np.array(
                [[1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0, 0.0],
                 [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  1.0, 0.0, 0.0, 0.0],
                 [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 1.0, 1.0],
                 [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                  0.0, 0.0, 1.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 1.0, 1.0]]
            )
        else:
            possible_rfs = np.array(np.zeros((6, 16)))
        ix = int((tstep / 1200) % 6)
        return possible_rfs[ix] * self.rand_mag

    # -----------------------------------------------------------------
    def set_next_state(
            self,
            next_state_array,
            close_point_smoothness_factors,
            x_cils,
            x_coas,
    ):
        self.curr_tstep += 1
        self.insert_state_array_into_system_history(next_state_array)

        verts = self.curr_state[:, [parameters.x_ix, parameters.y_ix]]

        edge_displacement_vectors_plus = geometry.calculate_edge_vectors(
            verts)
        edge_lengths = geometry.calculate_2D_vector_mags(
            edge_displacement_vectors_plus)

        self.curr_state[:, parameters.edge_lengths_ix] = edge_lengths
        self.curr_state[:, parameters.old_x_coa_ix] = x_coas
        self.curr_state[:, parameters.old_x_cil_ix] = x_cils

        # ==================================
        if self.curr_tstep == self.next_randomization_event_tpoint:
            self.next_randomization_event_tpoint += 1200

            # randomization event has occurred, so renew Rac kgtp rate
            # multipliers
            self.rac_rands = \
                self.renew_rac_rands(
                    self.curr_tstep)

        # store the Rac randomization factors for this timestep
        self.curr_state[:, parameters.rac_rands_ix] = self.rac_rands

        # ==================================
        rac_acts = self.curr_state[:, parameters.rac_act_ix]
        rho_acts = self.curr_state[:, parameters.rho_act_ix]

        self.curr_state[:, parameters.old_x_cil_ix] = x_cils
        self.curr_state[:, parameters.old_x_coa_ix] = x_coas

        rac_cyto = (
                1
                - np.sum(rac_acts)
                - np.sum(
            self.curr_state[:, parameters.rac_inact_ix]))
        rho_cyto = (
                1
                - np.sum(rho_acts)
                - np.sum(
            self.curr_state[:, parameters.rho_inact_ix]))

        insertion_array = np.zeros(16)
        insertion_array[0] = 1

        self.curr_state[:, parameters.rac_cyto_ix] = (
                rac_cyto * insertion_array)
        self.curr_state[:, parameters.rho_cyto_ix] = (rho_cyto *
                                                      insertion_array)

        sum_forces, edge_forces_plus, edge_forces_minus, uevs, rgtp_forces, \
        cyto_forces, \
        edge_strains, local_strains, unit_inside_pointing_vecs = \
            mechanics.calculate_forces(
                verts,
                rac_acts,
                rho_acts,
                self.rest_edge_len,
                self.stiffness_edge,
                self.halfmax_vertex_rgtp_act,
                self.const_protrusive,
                self.const_retractive,
                self.rest_area,
                self.stiffness_cyto,
            )

        self.curr_state[:, parameters.local_strains_ix] = local_strains

        edgeplus_lengths = geometry.calculate_edgeplus_lengths(verts)
        avg_edge_lengths = geometry.calculate_average_edge_length_around_nodes(
            edgeplus_lengths
        )

        conc_rac_acts = chemistry.calc_concs(
            rac_acts, avg_edge_lengths
        )

        conc_rho_acts = chemistry.calc_concs(
            rho_acts, avg_edge_lengths
        )

        self.curr_state[:, parameters.rac_rands_ix] = self.rac_rands

        kgtp_racs = chemistry.calculate_kgtp_rac(
            conc_rac_acts,
            self.halfmax_vertex_rgtp_conc,
            self.kgtp_rac,
            self.kgtp_rac_auto,
            x_coas,
            self.rac_rands,
            x_cils,
            close_point_smoothness_factors,
        )

        kdgtps_rac = chemistry.calculate_kdgtp_rac(
            conc_rho_acts,
            self.halfmax_vertex_rgtp_conc,
            self.kdgtp_rac,
            self.kdgtp_rho_on_rac,
            x_cils,
            self.halfmax_tension_inhib,
            self.tension_inhib,
            np.array([ls if ls > 0 else 0.0 for ls in local_strains]),
        )

        kdgtps_rho = chemistry.calculate_kdgtp_rho(
            conc_rac_acts,
            self.halfmax_vertex_rgtp_conc,
            self.kdgtp_rho,
            self.kdgtp_rac_on_rho,
        )

        kgtp_rhos = chemistry.calculate_kgtp_rho(
            conc_rho_acts,
            x_cils,
            self.halfmax_vertex_rgtp_conc,
            self.kgtp_rho,
            self.kgtp_rho_auto,
        )

        self.curr_state[:, parameters.kgtp_rac_ix] = kgtp_racs
        self.curr_state[:, parameters.kgtp_rho_ix] = kgtp_rhos

        self.curr_state[:, parameters.kdgtp_rac_ix] = kdgtps_rac

        self.curr_state[:, parameters.kdgtp_rho_ix] = kdgtps_rho

        # update mechanics parameters
        self.curr_state[:, [parameters.sum_forces_x_ix,
                            parameters.sum_forces_y_ix]] = sum_forces
        self.curr_state[:, [parameters.edge_forces_plus_x_ix,
                            parameters.edge_forces_plus_y_ix]] = \
            edge_forces_plus
        self.curr_state[:, [parameters.edge_forces_minus_x_ix,
                            parameters.edge_forces_minus_y_ix]] = \
            edge_forces_minus
        self.curr_state[:, [parameters.rgtp_forces_x_ix,
                            parameters.rgtp_forces_y_ix]] = rgtp_forces
        self.curr_state[:, [parameters.cyto_forces_x_ix,
                            parameters.cyto_forces_y_ix]] = cyto_forces

        self.curr_state[:, [parameters.uiv_x_ix,
                            parameters.uiv_y_ix]] = unit_inside_pointing_vecs
        self.curr_verts = verts

    # -----------------------------------------------------------------
    def pack_rhs_arguments(
            self,
            all_cells_verts,
            close_point_smoothness_factors,
            x_cils, x_coas,
    ):

        # print("-----------------------------------")
        # np.set_printoptions(suppress=True)
        # print("tstep: {}, cell: {}".format(self.curr_tstep, self.cell_ix))
        coa_update = np.abs(self.curr_state[:, parameters.old_x_coa_ix] -
                            x_coas) > 1e-4
        # if np.any(coa_update):
        #     print("old_coa = {}".format([float(x) if not x.is_integer() else
        #                                  int(x) for x in
        #                                  np.round(self.curr_state[:,
        #                                           parameters.old_x_coa_ix], 4)]))
        #     print("new_coa = {}".format([float(x) if not x.is_integer() else
        #                                  int(x) for x in np.round(x_coas, 4)]))
        # else:
        #     print("old_coa = no change")
        #     print("new_coa = no change")

        cil_update = np.abs(x_cils -
                            self.curr_state[:, parameters.old_x_cil_ix]) > 1e-4
        # if np.any(cil_update):
        #     print("old_cil = {}".format([float(x) if not x.is_integer() else
        #                                  int(x) for x in
        #                                  np.round(self.curr_state[:,
        #                                           parameters.old_x_cil_ix], 4)]))
        #     print("new_cil = {}".format([float(x) if not x.is_integer() else
        #                                  int(x) for x in np.round(x_cils, 4)]))
        # else:
        #     print("old_cil = no change")
        #     print("new_cil = no change")
        # print("-----------------------------------")

        return (
            self.num_cells,
            all_cells_verts,
            self.num_phase_vars,
            self.rac_act_ix,
            self.rest_edge_len,
            self.rac_inact_ix,
            self.rho_act_ix,
            self.rho_inact_ix,
            self.x_ix,
            self.y_ix,
            self.kgtp_rac,
            self.kdgtp_rac,
            self.kgtp_rho,
            self.kdgtp_rho,
            self.kgtp_rac_auto,
            self.kgtp_rho_auto,
            self.kdgtp_rho_on_rac,
            self.kdgtp_rac_on_rho,
            self.k_mem_on_vertex,
            self.k_mem_off,
            self.halfmax_vertex_rgtp_conc,
            self.diffusion_rgtp,
            self.vertex_eta,
            self.stiffness_edge,
            self.halfmax_vertex_rgtp_act,
            self.const_protrusive,
            self.const_retractive,
            self.rest_area,
            self.stiffness_cyto,
            x_coas,
            close_point_smoothness_factors,
            x_cils,
            self.halfmax_tension_inhib,
            self.tension_inhib,
            self.rac_rands,
            coa_update,
            cil_update,
        )

    # -----------------------------------------------------------------
    def execute_step(
            self,
            this_cell_ix,
            all_cells_verts,
            close_point_smoothness_factors,
            x_cils,
            x_coas,
            writer,
    ):
        dynamics.print_var = True

        num_cells = all_cells_verts.shape[0]

        state_array = dynamics.pack_state_array_from_system_history(
            self.phase_var_indices,
            self.ode_cellwide_phase_var_indices,
            self.curr_state,
        )

        rhs_args = self.pack_rhs_arguments(all_cells_verts,
                                           close_point_smoothness_factors,
                                           x_cils, x_coas)

        output_array = dynamics.eulerint(
            dynamics.cell_dynamics, state_array, np.array([0, 1]),
            rhs_args, self.num_int_steps, self.cell_ix, self.curr_tstep,
            self.rac_act_ix, self.rac_inact_ix,
            self.rho_act_ix, self.rho_inact_ix,
            self.x_ix, self.y_ix, writer)

        next_state_array = output_array[1]

        self.set_next_state(
            next_state_array,
            close_point_smoothness_factors,
            x_cils,
            x_coas,
        )

# ===============================================================
