import copy
import time

import numba as nb
import numpy as np

import cell
import geometry
import parameters

import hardio
from hardio import Writer
import chemistry

"""
Environment of cells.s
"""

MODE_EXECUTE = 0
MODE_OBSERVE = 1


@nb.jit(nopython=True)
def custom_floor(fp_number, roundoff_distance):
    a = int(fp_number)
    b = a + 1

    if abs(fp_number - b) < roundoff_distance:
        return b
    else:
        return a


# -----------------------------------------------------------------


def calc_bb_centre(bbox):
    x = bbox[0] + (bbox[1] - bbox[0]) * 0.5
    y = bbox[2] + (bbox[3] - bbox[2]) * 0.5

    return np.array([x, y])


# -----------------------------------------------------------------


class Environment:
    """Implementation of coupled map lattice model of a cell.
    """

    def __init__(
            self, out_dir, name, log_level, params
    ):
        self.params = params

        self.curr_tpoint = 0
        self.t = params["t"]
        self.dt = 1.0
        self.num_tsteps = int(np.ceil(params["final_t"]/self.t))
        self.num_tpoints = self.num_tsteps + 1
        self.timepoints = np.arange(0, self.num_tpoints)
        self.num_int_steps = params["num_int_steps"]
        self.num_cells = params["num_cells"]
        self.snap_period = params["snap_period"]

        self.writer = Writer(out_dir, name, log_level, params, self.snap_period)
        self.writer.save_header(params)

        self.env_cells = self.make_cells()
        self.cell_indices = np.arange(self.num_cells)
        self.exec_orders = np.zeros(
            (self.num_tpoints, self.num_cells), dtype=np.int64
        )

        self.all_geometry_tasks = np.array(
            geometry.create_dist_and_line_segment_intersection_test_args(
                self.num_cells, 16
            ),
            dtype=np.int64,
        )
        self.geometry_tasks_per_cell = np.array(
            [
                geometry.create_dist_and_line_segment_intersection_test_args_relative_to_specific_cell(
                    ci, self.num_cells, 16
                )
                for ci in range(self.num_cells)
            ],
            dtype=np.int64,
        )

        env_cells_verts = np.array(
            [x.curr_verts for x in self.env_cells]
        )
        cells_bb_array = \
            geometry.calc_init_cell_bbs(
                self.num_cells, env_cells_verts)
        (
            cells_node_distance_matrix,
            cells_line_segment_intersection_matrix,
        ) = geometry.init_lseg_intersects_and_dist_sq_matrices_old(
            self.num_cells,
            16,
            cells_bb_array,
            env_cells_verts,
        )
        self.x_cils_per_cell = []
        self.x_coas_per_cell = []
        self.update_coa_cil(env_cells_verts,
                            cells_node_distance_matrix,
                            cells_line_segment_intersection_matrix)
        for (ci, cell) in enumerate(self.env_cells):
            are_nodes_inside_other_cells = \
                geometry.check_if_nodes_inside_other_cells(
                    cell.cell_ix, 16, self.num_cells,
                    env_cells_verts)

            close_point_on_other_cells_to_each_node_exists, \
            close_point_on_other_cells_to_each_node, \
            close_point_on_other_cells_to_each_node_indices, \
            close_point_on_other_cells_to_each_node_projection_factors, \
            close_point_smoothness_factors = \
                geometry.do_close_points_to_each_node_on_other_cells_exist(
                    cell.num_cells, 16, cell.cell_ix,
                    env_cells_verts[cell.cell_ix],
                    cells_node_distance_matrix[ci],
                    cell.close_zero_at,
                    cell.close_one_at,
                    env_cells_verts, are_nodes_inside_other_cells)

            cell.initialize_cell(close_point_smoothness_factors,
                                 self.x_cils_per_cell[ci],
                                 self.x_coas_per_cell[ci])

        self.mode = MODE_EXECUTE
        self.animation_settings = None

    # -----------------------------------------------------------------

    def make_cells(self):
        env_cells, init_cell_bbs = self.create_cell_group(self.params)
        return np.array(env_cells)

    # -----------------------------------------------------------------

    def create_cell_group(
            self, params
    ):
        num_cells = params["num_cells"]
        cell_group_bb = params["cell_group_bbox"]
        cell_r = params["cell_r"]

        init_cell_bbs = self.calculate_cell_bbs(
            num_cells,
            cell_r,
            cell_group_bb,
        )

        cells_in_group = []

        for ci, bbox in enumerate(init_cell_bbs):
            (
                init_verts,
                rest_edge_len,
                rest_area,
            ) = self.create_default_init_cell_verts(
                bbox, cell_r
            )

            params.update(
                [
                    ("init_verts", init_verts),
                    ("rest_edge_len", rest_edge_len),
                    ("rest_area", rest_area),
                ]
            )

            new_cell = cell.Cell(
                0,
                ci,
                params,
                self.dt,
            )

            cells_in_group.append(new_cell)

        return cells_in_group, init_cell_bbs

    # -----------------------------------------------------------------

    @staticmethod
    def calculate_cell_bbs(
            num_cells,
            cell_r,
            cell_group_bb,
    ):

        cell_bbs = np.zeros((num_cells, 4), dtype=np.float64)
        xmin, xmax, ymin, ymax = cell_group_bb
        x_length = xmax - xmin
        y_length = ymax - ymin

        cell_diameter = 2 * cell_r

        # check if cells can fit in given bounding box
        total_cell_group_area = num_cells * (np.pi * cell_r ** 2)
        cell_group_bb_area = abs(x_length * y_length)

        if total_cell_group_area > cell_group_bb_area:
            raise Exception(
                "Cell group bounding box is not big enough to contain all "
                "cells given cell_r constraint."
            )
        num_cells_along_x = custom_floor(x_length / cell_diameter, 1e-6)
        num_cells_along_y = custom_floor(y_length / cell_diameter, 1e-6)

        cell_x_coords = (
                xmin +
                np.sign(x_length) *
                np.arange(num_cells_along_x) *
                cell_diameter)
        cell_y_coords = (
                ymin +
                np.sign(y_length) *
                np.arange(num_cells_along_y) *
                cell_diameter)
        x_step = np.sign(x_length) * cell_diameter
        y_step = np.sign(y_length) * cell_diameter

        xi = 0
        yi = 0
        for ci in range(num_cells):
            cell_bbs[ci] = [
                cell_x_coords[xi],
                cell_x_coords[xi] + x_step,
                cell_y_coords[yi],
                cell_y_coords[yi] + y_step,
            ]

            if yi == (num_cells_along_y - 1):
                yi = 0
                xi += 1
            else:
                yi += 1

        return cell_bbs

    # -----------------------------------------------------------------

    @staticmethod
    def create_default_init_cell_verts(
            bbox, cell_r
    ):
        cell_centre = calc_bb_centre(bbox)

        cell_node_thetas = np.pi * \
                           np.linspace(0, 2, endpoint=False, num=16)
        cell_verts = np.transpose(
            np.array(
                [
                    cell_r * np.cos(cell_node_thetas),
                    cell_r * np.sin(cell_node_thetas),
                ]
            )
        )

        # rotation_theta = np.random.rand()*2*np.pi
        # cell_verts = np.array([
        # geometry.rotate_2D_vector_CCW_by_theta(rotation_theta, x) for x in
        # cell_verts], dtype=np.float64)
        cell_verts = np.array(
            [[x + cell_centre[0], y + cell_centre[1]] for x, y in
             cell_verts],
            dtype=np.float64,
        )

        edge_vectors = geometry.calculate_edge_vectors(cell_verts)

        edge_lengths = geometry.calculate_vec_mags(edge_vectors)

        rest_edge_len = np.average(edge_lengths)

        rest_area = geometry.calculate_polygon_area(cell_verts)
        if rest_area < 0:
            raise Exception("Resting area was calculated to be negative.")
        return cell_verts, rest_edge_len, rest_area

    # -----------------------------------------------------------------
    def execute_system_dynamics_in_random_sequence(
            self,
            t,
            cells_node_distance_matrix,
            cells_bb_array,
            cells_line_segment_intersection_matrix,
            env_cells_verts,
            environment_cells,
    ):
        execution_sequence = self.cell_indices
        # np.random.shuffle(execution_sequence)

        self.exec_orders[t] = np.copy(execution_sequence)

        for ci in execution_sequence:
            current_cell = environment_cells[ci]
            x_cils = self.x_cils_per_cell[ci]
            x_coas = self.x_coas_per_cell[ci]

            are_nodes_inside_other_cells = \
                geometry.check_if_nodes_inside_other_cells(
                    current_cell.cell_ix, 16, self.num_cells,
                    env_cells_verts)

            close_point_on_other_cells_to_each_node_exists, \
            close_point_on_other_cells_to_each_node, \
            close_point_on_other_cells_to_each_node_indices, \
            close_point_on_other_cells_to_each_node_projection_factors, \
            close_point_smoothness_factors = \
                geometry.do_close_points_to_each_node_on_other_cells_exist(
                    current_cell.num_cells, 16, current_cell.cell_ix,
                    env_cells_verts[current_cell.cell_ix],
                    cells_node_distance_matrix[ci],
                    current_cell.close_zero_at,
                    current_cell.close_one_at,
                    env_cells_verts, are_nodes_inside_other_cells)

            self.writer.begin_cell_save()
            current_cell.execute_step(
                env_cells_verts,
                close_point_smoothness_factors,
                x_cils,
                x_coas,
                self.writer,
            )
            self.writer.finish_cell_save()

            this_cell_coords = current_cell.curr_verts

            env_cells_verts[ci] = this_cell_coords

            cells_bb_array[ci] = \
                geometry.calculate_polygon_bb(this_cell_coords)
            geometry.update_line_segment_intersection_and_dist_squared_matrices(
                4,
                self.geometry_tasks_per_cell[ci],
                env_cells_verts,
                cells_bb_array,
                cells_node_distance_matrix,
                cells_line_segment_intersection_matrix,
                sequential=True,
            )

            self.update_coa_cil(env_cells_verts, cells_node_distance_matrix,
                                cells_line_segment_intersection_matrix)

        return (
            cells_node_distance_matrix,
            cells_bb_array,
            cells_line_segment_intersection_matrix,
            env_cells_verts,
        )

    # -----------------------------------------------------------------

    def update_coa_cil(self, env_cells_verts, cells_node_distance_matrix,
                       cells_line_segment_intersection_matrix):
        self.x_cils_per_cell = []
        self.x_coas_per_cell = []
        for ci in range(self.num_cells):
            cell = self.env_cells[ci]
            x_coas = chemistry.calculate_x_coas(
                self.num_cells,
                ci,
                cell.coa_distrib_exp,
                cell.coa_mag,
                cells_node_distance_matrix[ci],
                cells_line_segment_intersection_matrix[ci],
                cell.coa_los_penalty,
            )

            are_nodes_inside_other_cells = \
                geometry.check_if_nodes_inside_other_cells(
                    cell.cell_ix, 16, self.num_cells,
                    env_cells_verts)

            close_point_on_other_cells_to_each_node_exists, \
            close_point_on_other_cells_to_each_node, \
            close_point_on_other_cells_to_each_node_indices, \
            close_point_on_other_cells_to_each_node_projection_factors, \
            close_point_smoothness_factors = \
                geometry.do_close_points_to_each_node_on_other_cells_exist(
                    cell.num_cells, 16, cell.cell_ix,
                    env_cells_verts[cell.cell_ix],
                    cells_node_distance_matrix[ci],
                    cell.close_zero_at,
                    cell.close_one_at,
                    env_cells_verts, are_nodes_inside_other_cells)
            x_cils = chemistry.calc_x_cils(cell.cell_ix, cell.num_cells,
                                           cell.cil_mag,
                                           close_point_smoothness_factors)
            self.x_coas_per_cell.append(x_coas)
            self.x_cils_per_cell.append(x_cils)

    def execute_system_dynamics(
            self,
    ):
        simulation_st = time.time()
        num_cells = self.num_cells

        environment_cells = self.env_cells
        all_cell_verts = np.array([x.curr_verts for x in environment_cells])

        cell_bbs = \
            geometry.calc_init_cell_bbs(
                num_cells, all_cell_verts)
        (
            cells_node_distance_matrix,
            cells_line_segment_intersection_matrix,
        ) = geometry.init_lseg_intersects_and_dist_sq_matrices_old(
            num_cells,
            16,
            cell_bbs,
            all_cell_verts,
        )

        cell_group_indices = []

        for a_cell in self.env_cells:
            cell_group_indices.append(a_cell.cell_group_ix)

        # self.writer.begin_initial_save()
        # for (ci, cell) in enumerate(self.env_cells):
        #     x_cils = self.x_cils_per_cell[ci]
        #     x_coas = self.x_coas_per_cell[ci]
        #     cell.save_with_writer(x_cils, x_coas, self.writer)
        # self.writer.finish_initial_save()

        for t in self.timepoints:
            self.writer.begin_tpoint_save(t)
            (
                cells_node_distance_matrix,
                cell_bbs,
                cells_line_segment_intersection_matrix,
                all_cell_verts,
            ) = self.execute_system_dynamics_in_random_sequence(
                t,
                cells_node_distance_matrix,
                cell_bbs,
                cells_line_segment_intersection_matrix,
                all_cell_verts,
                environment_cells,
            )
            self.writer.finish_tpoint_save()

        simulation_et = time.time()
        self.writer.finish()

        simulation_time = np.round(simulation_et - simulation_st, decimals=2)

        print(
            ("Time taken to complete simulation: {}s".format(simulation_time))
        )

# -----------------------------------------------------------------
