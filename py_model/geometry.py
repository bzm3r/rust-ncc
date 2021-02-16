# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:22:06 2015

@author: Brian
"""

import numpy as np
import math
import numba as nb
import threading
from hardio import GEOM_ALG_INFO


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def rotate_vector_ccw(vector):
    x, y = vector

    result_vector = np.empty(2, dtype=np.float64)

    result_vector[0] = -1.0 * y
    result_vector[1] = x

    return result_vector


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def rotate_vectors_ccw(vectors):
    num_vectors = vectors.shape[0]
    rotated_vectors = np.empty_like(vectors)

    for i in range(num_vectors):
        rotated_vector = rotate_vector_ccw(vectors[i])
        for j in range(2):
            rotated_vectors[i, j] = rotated_vector[j]

    return rotated_vectors


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_vec_mag(vector):
    x = vector[0]
    y = vector[1]

    return math.sqrt(x * x + y * y)


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_vec_mags(vectors):
    num_vectors = vectors.shape[0]
    vector_mags = np.empty(num_vectors)

    for i in range(num_vectors):
        vector_mags[i] = calculate_vec_mag(vectors[i])

    return vector_mags


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_vec_directions(vector):
    x = vector[0]
    y = vector[1]

    abs_x = abs(x)
    abs_y = abs(y)

    if abs_x < 1e-4:
        if abs_y < 1e-4:
            return -1.0
        else:
            phi = np.pi / 2
    else:
        phi = np.arctan(abs_y / abs_x)

    if x >= 0:
        if y >= 0:
            return phi
        else:
            return 2 * np.pi - phi
    else:
        if y >= 0:
            return np.pi - phi
        else:
            return np.pi + phi


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calc_vec_mag(vector):
    normalized_vector = np.empty(2, dtype=np.float64)

    mag = calculate_vec_mag(vector)
    x, y = vector

    if mag < 1e-4:
        normalized_vector[0] = np.nan
        normalized_vector[1] = np.nan
    else:
        normalized_vector[0] = x / mag
        normalized_vector[1] = y / mag

    return normalized_vector


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def normalize_vectors(vectors):
    num_vectors = vectors.shape[0]
    normalized_vectors = np.empty((num_vectors, 2), dtype=np.float64)

    for i in range(num_vectors):
        vector = vectors[i]

        normalized_vector = calc_vec_mag(vector)

        if normalized_vector[0] == np.nan or normalized_vector[1] == np.nan:
            normalized_vectors[i, 0] = np.nan
            normalized_vectors[i, 1] = np.nan
        else:
            normalized_vectors[i, 0] = normalized_vector[0]
        normalized_vectors[i, 1] = normalized_vector[1]

    return normalized_vectors


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_edge_vectors(verts):
    num_vertices = verts.shape[0]
    edge_displacement_vectors = np.empty((num_vertices, 2), dtype=np.float64)
    for i in range(num_vertices):
        i_plus1 = (i + 1) % num_vertices

        x_a, y_a = verts[i]
        x_b, y_b = verts[i_plus1]

        edge_displacement_vectors[i, 0] = x_b - x_a
        edge_displacement_vectors[i, 1] = y_b - y_a

    return edge_displacement_vectors


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def get_vector_normals(num_vectors, vectors):
    normal_to_vectors = np.empty((num_vectors, 2), dtype=np.float64)

    for i in range(num_vectors):
        x, y = vectors[i]
        normal_to_vectors[i, 0] = -1 * y
        normal_to_vectors[i, 1] = x

    return normal_to_vectors


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_unit_inside_pointing_vecs(verts):
    num_vertices = verts.shape[0]

    unit_inside_pointing_vecs = np.empty((num_vertices, 2), dtype=np.float64)

    tangent_vector = np.empty(2, dtype=np.float64)

    for i in range(num_vertices):
        i_plus_1 = (i + 1) % num_vertices

        edge_vector_to_plus = verts[i_plus_1] - verts[i]
        edge_vector_to_plus_normalized = calc_vec_mag(
            edge_vector_to_plus)

        i_minus_1 = (i - 1) % num_vertices

        edge_vector_from_minus = verts[i] - verts[i_minus_1]
        edge_vectors_from_minus_normalized = calc_vec_mag(
            edge_vector_from_minus)

        tangent_vector[0] = (
                edge_vector_to_plus_normalized[0] +
                edge_vectors_from_minus_normalized[0]
        )
        tangent_vector[1] = (
                edge_vector_to_plus_normalized[1] +
                edge_vectors_from_minus_normalized[1]
        )

        tangent_vector_normalized = calc_vec_mag(tangent_vector)

        x_part, y_part = rotate_vector_ccw(tangent_vector_normalized)
        unit_inside_pointing_vecs[i, 0] = x_part
        unit_inside_pointing_vecs[i, 1] = y_part

    return unit_inside_pointing_vecs


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_unit_inside_pointing_vecs_per_timestep(verts_per_timestep):
    num_tsteps = verts_per_timestep.shape[0]
    num_vertices = verts_per_timestep.shape[1]

    unit_inside_pointing_vecs_per_timestep = np.empty(
        (num_tsteps, num_vertices, 2), dtype=np.float64
    )

    for t in range(num_tsteps):
        unit_inside_pointing_vecs_per_timestep[
            t] = calculate_unit_inside_pointing_vecs(
            verts_per_timestep[t]
        )

    return unit_inside_pointing_vecs_per_timestep


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_tangent_vecs_and_inside_pointing_vecs(num_vertices, verts):
    unit_inside_pointing_vecs = np.empty((num_vertices, 2), dtype=np.float64)
    tangent_vectors = np.empty((num_vertices, 2), dtype=np.float64)

    tangent_vector = np.empty(2, dtype=np.float64)
    for i in range(num_vertices):
        i_plus_1 = (i + 1) % num_vertices

        edge_vector_to_plus = calculate_vector_from_p1_to_p2_given_vectors(
            verts[i], verts[i_plus_1]
        )
        edge_vector_to_plus_normalized = calc_vec_mag(
            edge_vector_to_plus)

        i_minus_1 = (i - 1) % num_vertices

        edge_vector_from_minus = calculate_vector_from_p1_to_p2_given_vectors(
            verts[i_minus_1], verts[i]
        )
        edge_vectors_from_minus_normalized = calc_vec_mag(
            edge_vector_from_minus)

        tangent_vector[0] = (
                edge_vector_to_plus_normalized[0] +
                edge_vectors_from_minus_normalized[0]
        )
        tangent_vector[1] = (
                edge_vector_to_plus_normalized[1] +
                edge_vectors_from_minus_normalized[1]
        )

        tangent_vector_normalized = calc_vec_mag(tangent_vector)
        tangent_vectors[i, 0], tangent_vectors[i, 1] = tangent_vector_normalized

        normal_to_tangent_vector = rotate_vector_ccw(
            tangent_vector_normalized)
        unit_inside_pointing_vecs[i, 0], unit_inside_pointing_vecs[
            i, 1
        ] = normal_to_tangent_vector

    return unit_inside_pointing_vecs, tangent_vectors


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def check_sign_projection_of_a_on_b(a, b):
    """
    The magnitude of the projection of a on b is given by:
        (a \\dot b)/|b|

    However, |b| is always >= 0.

    Thus, the sign is determined solely by the sign of a \\dot b.
    """
    if (a[0] * b[0] + a[1] * b[1]) > 0:
        return 1
    else:
        return 0


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def is_left(p0, p1, p2):
    """
    Input:  three points P0, P1, and P2
    Return: > 0 for P2 left of the line through P0 to P1
            = 0 for P2 on the line
            < 0 for P2 right of the line
    """
    p0x, p0y = p0
    p1x, p1y = p1
    p2x, p2y = p2

    # logging.log(level=GEOM_ALG_INFO,
    #             msg=(p1x - p0x) * (p2y - p0y) - (p2x - p0x) * (p1y - p0y))
    return (p1x - p0x) * (p2y - p0y) - (p2x - p0x) * (p1y - p0y)


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_polygon_bb(polygon):
    num_vertices = polygon.shape[0]

    min_x = 0.0
    max_x = 0.0

    min_y = 0.0
    max_y = 0.0

    for i in range(num_vertices):
        this_x, this_y = polygon[i]
        if i == 0:
            min_x = this_x
            max_x = min_x

            min_y = this_y
            max_y = min_y
        else:
            if this_x < min_x:
                min_x = this_x
            elif this_x > max_x:
                max_x = this_x

            if this_y < min_y:
                min_y = this_y
            elif this_y > max_y:
                max_y = this_y

    return min_x, max_x, min_y, max_y


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_bb_polygon(polygon):
    num_vertices = polygon.shape[0]

    min_x = 0.0
    max_x = 0.0

    min_y = 0.0
    max_y = 0.0

    for i in range(num_vertices):
        if i == 0:
            min_x = polygon[i][0]
            max_x = min_x

            min_y = polygon[i][1]
            max_y = min_y
        else:
            this_x = polygon[i][0]
            this_y = polygon[i][1]

            if this_x < min_x:
                min_x = this_x
            elif this_x > max_x:
                max_x = this_x

            if this_y < min_y:
                min_y = this_y
            elif this_y > max_y:
                max_y = this_y

    bb_polygon = np.zeros((4, 2), dtype=np.float64)

    for i in range(4):
        if i < 2:
            bb_polygon[i][1] = min_y
            if i % 2 == 0:
                bb_polygon[i][0] = min_x
            else:
                bb_polygon[i][0] = max_x
        else:
            bb_polygon[i][1] = max_y
            if i % 2 == 1:
                bb_polygon[i][0] = min_x
            else:
                bb_polygon[i][0] = max_x

    return bb_polygon


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def calc_init_cell_bbs(
        num_cells, environment_cells_verts
):
    bb_polygon_array = np.zeros((num_cells, 4), dtype=np.float64)

    for ci in range(num_cells):
        bb_polygon_array[ci] = calculate_polygon_bb(
            environment_cells_verts[ci]
        )

    return bb_polygon_array


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def is_point_in_polygon_bb_given_bb(
        test_point, min_x, max_x, min_y, max_y
):
    tp_x = test_point[0]
    tp_y = test_point[1]

    return (min_x < tp_x < max_x) and (min_y < tp_y < max_y)


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def is_point_in_polygon_bb(test_point, polygon):
    min_x, max_x, min_y, max_y = calculate_polygon_bb(polygon)

    return is_point_in_polygon_bb_given_bb(
        test_point, min_x, max_x, min_y, max_y
    )


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def is_point_in_polygon_without_bb_check(test_point, polygon):
    num_vertices = polygon.shape[0]
    wn = 0
    test_point_y = np.around(test_point[1], 4)

    # count number of intersections of positive-x direction ray emanating from test_point with polygon edges
    for i in range(num_vertices):
        p_start = polygon[i]
        p_end = polygon[(i + 1) % num_vertices]

        p_start_y = np.around(p_start[1], 4)

        p_end_y = np.around(p_end[1], 4)

        # logging.log(level=GEOM_ALG_INFO, msg="abs(p_start.y - test_point.y): {}".format(
        #     np.abs(p_start_y -
        #            test_point_y)))
        # logging.log(level=GEOM_ALG_INFO,
        #             msg="psy: {}, pey: {}, tpy: {}".format(p_start_y, p_end_y,
        #                                                    test_point_y))
        psy_eq_tpy = np.abs(p_start_y - test_point_y) < 1e-8

        if psy_eq_tpy and (np.abs(p_start[0] - test_point[0]) < 1e-8):
            continue
        psy_leq_tpy = psy_eq_tpy or p_start_y < test_point_y  # p_start_y <= test_point_y
        tpy_le_pey = test_point_y < p_end_y
        cond0 = psy_leq_tpy and tpy_le_pey

        # logging.log(level=GEOM_ALG_INFO, msg="psy == tpy: {}".format(psy_eq_tpy))
        # logging.log(level=GEOM_ALG_INFO, msg="psy <= tpy: {}".format(psy_leq_tpy))
        # logging.log(level=GEOM_ALG_INFO, msg="pey < tpy: {}".format(tpy_le_pey))
        # logging.log(level=GEOM_ALG_INFO, msg="psy <= tpy < pey: {}".format(cond0))

        pey_le_tpy = p_end_y < test_point_y
        tpy_leq_psy = psy_eq_tpy or test_point_y < p_start_y  # test_point_y <= p_start_y
        cond1 = pey_le_tpy and tpy_leq_psy

        # logging.log(level=GEOM_ALG_INFO, msg="psy == tpy: {}".format(psy_eq_tpy))
        # logging.log(level=GEOM_ALG_INFO, msg="tpy <= py: {}".format(tpy_leq_psy))
        # logging.log(level=GEOM_ALG_INFO, msg="pey < tpy: {}".format(pey_le_tpy))
        # logging.log(level=GEOM_ALG_INFO, msg="pey < tpy <= psy: {}".format(cond1))
        #
        # logging.log(level=97,
        #             msg="ovi: {}, cond0: {}, cond1: {}".format(i, cond0, cond1))
        if cond0:
            # upward crossing
            is_tp_left_of_edge = is_left(p_start, p_end, test_point)

            if is_tp_left_of_edge > 0:
                # logging.log(level=GEOM_ALG_INFO, msg="upward crossing")
                # positive-x direction ray emanating from test_point wil intersect with this edge if left of it
                wn += 1
        elif cond1:
            # downward crossing
            is_tp_left_of_edge = is_left(p_start, p_end, test_point)

            if is_tp_left_of_edge < 0:
                # logging.log(level=GEOM_ALG_INFO, msg="downward crossing")
                # positive-x direction ray emanating from test_point will intersect with this edge if left of it
                wn -= 1

    # logging.log(level=GEOM_ALG_INFO, msg="wn = {}, wn != 0 == {}".format(wn, wn != 0))
    return wn != 0


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def is_point_in_polygon(test_point, polygon):
    is_test_point_in_poly_bb = is_point_in_polygon_bb(test_point, polygon)

    if is_test_point_in_poly_bb == 0:
        return 0
    else:
        return is_point_in_polygon_without_bb_check(test_point, polygon)


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def is_point_in_polygon_given_bb(
        test_point, polygon, min_x, max_x, min_y, max_y
):
    is_test_point_in_poly_bb = is_point_in_polygon_bb_given_bb(
        test_point, min_x, max_x, min_y, max_y
    )

    return is_test_point_in_poly_bb and \
           is_point_in_polygon_without_bb_check(test_point, polygon)


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def is_point_in_polygons(point, polygons):
    num_polygons = polygons.shape[0]
    if num_polygons != 0:
        for pi in range(num_polygons):
            polygon = polygons[pi]
            if is_point_in_polygon(point, polygon.shape[0], polygon):
                return 1
            else:
                continue
    else:
        return 0

    return 0


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def are_points_inside_polygon(points, polygon):
    num_points = points.shape[0]
    results = np.zeros(num_points, dtype=np.bool)

    # min_x, max_x, min_y, max_y = calculate_polygon_bb(polygon)

    for index in range(num_points):
        test_point = points[index]
        # logging.log(level=97, msg="testing vertex: {}".format(index))
        results[index] = is_point_in_polygon_without_bb_check(test_point,
                                                              polygon)

    return results


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def are_points_inside_polygons(points, polygons):
    num_points = points.shape[0]
    results = np.zeros(num_points, dtype=np.int64)

    for index in range(num_points):
        test_point = points[index]
        for polygon in polygons:
            min_x, max_x, min_y, max_y = calculate_polygon_bb(polygon)
            if (
                    is_point_in_polygon_given_bb(
                        test_point, polygon, min_x, max_x, min_y, max_y
                    )
                    == 1
            ):
                results[index] = 1
                continue

    return results


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def check_if_nodes_inside_other_cells(
        this_cell_ix, num_vertices, num_cells, all_cells_verts
):
    are_nodes_inside_other_cells = np.zeros((num_vertices, num_cells),
                                            dtype=np.int64)

    this_cell_verts = all_cells_verts[this_cell_ix]

    for other_ci in range(num_cells):
        if other_ci != this_cell_ix:
            are_nodes_inside_current_cell = are_points_inside_polygon(
                this_cell_verts, all_cells_verts[other_ci]
            )
            for ni in range(num_vertices):
                are_nodes_inside_other_cells[ni][
                    other_ci
                ] = are_nodes_inside_current_cell[ni]

    return are_nodes_inside_other_cells


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_polygon_area(vertex_coords):
    """
    http://geomalgorithms.com/a01-_area.html
    """
    num_vertices = vertex_coords.shape[0]

    area = 0
    for i in range(num_vertices):
        j = (i + 1) % num_vertices
        k = (i - 1) % num_vertices
        area += vertex_coords[i, 0] * (
                vertex_coords[j, 1] - vertex_coords[k, 1])

    return area * 0.5


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_edgeplus_lengths(verts):
    num_vertices = verts.shape[0]
    edgeplus_lengths = np.zeros(num_vertices, dtype=np.float64)

    for ni in range(num_vertices):
        this_vert = verts[ni]
        plus_vert = verts[(ni + 1) % num_vertices]

        edgeplus_lengths[ni] = calculate_dist_between_points_given_vectors(
            this_vert, plus_vert
        )

    return edgeplus_lengths


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_average_edge_length_around_nodes(edgeplus_lengths):
    num_vertices = edgeplus_lengths.shape[0]
    avg_edge_lengths = np.zeros(num_vertices, dtype=np.float64)

    for ni in range(num_vertices):
        this_node_edgeplus_length = edgeplus_lengths[ni]
        last_node_edgeplus_length = edgeplus_lengths[(ni - 1) % num_vertices]

        avg_edge_lengths[ni] = (
                0.5 * this_node_edgeplus_length + 0.5 * last_node_edgeplus_length
        )

    return avg_edge_lengths


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_perimeter(verts):
    num_vertices = verts.shape[0]
    edgeplus_lengths = calculate_edgeplus_lengths(verts)

    perimeter = 0.0
    for ni in range(num_vertices):
        perimeter += edgeplus_lengths[ni]

    return perimeter


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def multiply_vectors_by_scalars(vectors, scalars):
    num_vectors = vectors.shape[0]
    result = np.empty((num_vectors, 2))
    for i in range(num_vectors):
        x, y = vectors[i]
        scalar = scalars[i]
        result[i, 0] = scalar * x
        result[i, 1] = scalar * y

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def multiply_vectors_by_scalar(vectors, scalar):
    num_vectors = vectors.shape[0]
    result = np.empty((num_vectors, 2))

    for i in range(num_vectors):
        x, y = vectors[i]
        result[i, 0] = scalar * x
        result[i, 1] = scalar * y

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_squared_dist(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]

    return dx * dx + dy * dy


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def create_initial_distance_squared_matrix(
        num_cells, num_vertices_per_cell, init_all_cells_verts
):
    distance_squared_matrix = -1 * np.ones(
        (num_cells, num_vertices_per_cell, num_cells, num_vertices_per_cell),
        dtype=np.float64,
    )
    info_update_tracker = np.zeros_like(distance_squared_matrix, dtype=np.int64)

    for ci in range(num_cells):
        for ni in range(num_vertices_per_cell):
            for other_ci in range(num_cells):
                if ci != other_ci:
                    relevant_info_update_tracker_slice = \
                        info_update_tracker[ci][ni][
                            other_ci
                        ]
                    for other_ni in range(num_vertices_per_cell):
                        if relevant_info_update_tracker_slice[other_ni] != 1:
                            info_update_tracker[ci][ni][other_ci][other_ni] = 1
                            info_update_tracker[other_ci][other_ni][ci][ni] = 1

                            this_node = init_all_cells_verts[ci][ni]
                            other_node = init_all_cells_verts[other_ci][
                                other_ni]

                            squared_dist = calculate_squared_dist(this_node,
                                                                  other_node)
                            distance_squared_matrix[ci][ni][other_ci][
                                other_ni
                            ] = squared_dist
                            distance_squared_matrix[other_ci][other_ni][ci][
                                ni
                            ] = squared_dist

    return distance_squared_matrix


# -----------------------------------------------------------------
@nb.jit(nopython=True, parallel=True)
def update_distance_squared_matrix_old(
        last_updated_cell_ix,
        num_cells,
        num_vertices_per_cell,
        all_cells_verts,
        distance_squared_matrix,
):
    # 1 if info has been updated, 0 otherwise
    info_update_tracker = np.zeros_like(distance_squared_matrix, dtype=np.int64)

    for ni in range(num_vertices_per_cell):
        for other_ci in range(num_cells):
            if other_ci != last_updated_cell_ix:
                relevant_info_update_tracker_slice = info_update_tracker[
                    last_updated_cell_ix
                ][ni][other_ci]
                for other_ni in range(num_vertices_per_cell):
                    if relevant_info_update_tracker_slice[other_ni] != 1:
                        info_update_tracker[last_updated_cell_ix][ni][other_ci][
                            other_ni
                        ] = 1
                        info_update_tracker[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = 1

                        this_node = all_cells_verts[last_updated_cell_ix][ni]
                        other_node = all_cells_verts[other_ci][other_ni]

                        squared_dist = calculate_squared_dist(this_node,
                                                              other_node)
                        distance_squared_matrix[last_updated_cell_ix][ni][
                            other_ci][
                            other_ni
                        ] = squared_dist
                        distance_squared_matrix[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = squared_dist

    return distance_squared_matrix


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_vector_from_p1_to_p2_given_vectors(p1, p2):
    p1x = p1[0]
    p1y = p1[1]

    p2x = p2[0]
    p2y = p2[1]

    result = np.empty(2, dtype=np.float64)

    result[0] = p2x - p1x
    result[1] = p2y - p1y

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_vector_from_p1_to_p2_given_coords(p1x, p1y, p2x, p2y):
    result = np.empty(2, dtype=np.float64)

    result[0] = p2x - p1x
    result[1] = p2y - p1y

    return result


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_projection_of_a_on_b(a, b):
    b_mag = calculate_vec_mag(b)

    return (a[0] * b[0] + a[1] * b[1]) / b_mag


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_dist_between_points_given_vectors(p1, p2):
    displacement_vector = p2 - p1
    distance = calculate_vec_mag(displacement_vector)

    return distance


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_dist_between_points_given_coords(p1x, p1y, p2x, p2y):
    displacement_vector = calculate_vector_from_p1_to_p2_given_coords(
        p1x, p1y, p2x, p2y
    )
    distance = calculate_vec_mag(displacement_vector)

    return distance


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def find_closest_node_on_other_cells_for_each_node_on_this_cell(
        num_cells, num_vertices, this_ci, dist_squared_array
):
    closest_nodes_on_other_cells = -1 * np.ones(
        (num_vertices, num_cells), dtype=np.int64
    )

    for ni in range(num_vertices):
        this_node_dist_squared_array = dist_squared_array[ni]
        for ci in range(num_cells):
            if ci != this_ci:
                closest_node = -1
                closest_node_dist = -1

                this_node_dist_squared_wrt_other_cell = \
                    this_node_dist_squared_array[ci]
                for other_ni in range(num_vertices):
                    dist_bw_nodes = this_node_dist_squared_wrt_other_cell[
                        other_ni]

                    if other_ni == 0 or dist_bw_nodes < closest_node_dist:
                        closest_node = other_ni
                        closest_node_dist = dist_bw_nodes

                closest_nodes_on_other_cells[ni][ci] = closest_node

    return closest_nodes_on_other_cells


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def calculate_closest_point_dist_squared(
        num_vertices, this_nc, other_cell_verts, closest_ni
):
    closest_point_node_indices = np.zeros(2, dtype=np.int64)
    closest_point_node_indices[0] = closest_ni
    closest_point_node_indices[1] = -1

    closest_ni_plus1 = (closest_ni + 1) % num_vertices
    closest_ni_minus1 = (closest_ni - 1) % num_vertices

    closest_nc = other_cell_verts[closest_ni]
    closest_nc_plus1 = other_cell_verts[closest_ni_plus1]
    closest_nc_minus1 = other_cell_verts[closest_ni_minus1]

    closest_to_this = this_nc - closest_nc
    plus1_vector = closest_nc_plus1 - closest_nc
    minus1_vector = closest_nc_minus1 - closest_nc

    proj_for_plus1 = calculate_projection_of_a_on_b(closest_to_this,
                                                    plus1_vector)
    if 0 < proj_for_plus1 < 1:
        closest_pc = closest_nc + proj_for_plus1 * plus1_vector
        closest_point_node_indices[1] = closest_ni_plus1
        return (
            calculate_squared_dist(this_nc, closest_pc),
            closest_pc,
            closest_point_node_indices,
            proj_for_plus1,
        )

    proj_for_minus1 = calculate_projection_of_a_on_b(closest_to_this,
                                                     minus1_vector)
    if 0 < proj_for_minus1 < 1:
        closest_pc = closest_nc + proj_for_minus1 * minus1_vector
        closest_point_node_indices[1] = closest_ni_minus1
        return (
            calculate_squared_dist(this_nc, closest_pc),
            closest_pc,
            closest_point_node_indices,
            proj_for_minus1,
        )

    return -1, np.zeros(2, dtype=np.float64), closest_point_node_indices, 1.0


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def closeness_smoothening_linear_function(zero_until, one_at, x):
    m = 1.0 / (one_at - zero_until)
    b = -zero_until * m

    if x > zero_until:
        return 0.0
    elif x < one_at:
        return 1.0
    else:
        return m * x + b


# @nb.jit(nopython=True)
def do_close_points_to_each_node_on_other_cells_exist(
        num_cells,
        num_vertices,
        this_ci,
        this_cell_verts,
        dist_squared_array,
        close_0_at,
        close_1_at,
        all_cells_verts,
        are_nodes_inside_other_cells,
):
    close_points_exist = np.zeros((num_vertices, num_cells), dtype=np.int64)
    close_points = np.zeros((num_vertices, num_cells, 2), dtype=np.float64)
    close_points_node_indices = np.zeros((num_vertices, num_cells, 2),
                                         dtype=np.int64)
    close_points_node_projection_factors = np.ones(
        (num_vertices, num_cells), dtype=np.float64
    )
    close_point_smoothness_factors = np.zeros(
        (num_vertices, num_cells), dtype=np.float64
    )

    closest_nodes_on_other_cells = find_closest_node_on_other_cells_for_each_node_on_this_cell(
        num_cells, num_vertices, this_ci, dist_squared_array
    )

    close_0 = np.sqrt(close_0_at)
    close_1 = np.sqrt(close_1_at)

    for ni in range(num_vertices):
        closest_nodes_to_this_node = closest_nodes_on_other_cells[ni]
        relevant_dist_squared_array_slice = dist_squared_array[ni]
        this_nc = this_cell_verts[ni]

        for ci in range(num_cells):
            if ci != this_ci:
                closest_ni = closest_nodes_to_this_node[ci]
                closest_node_dist = relevant_dist_squared_array_slice[ci][
                    closest_ni]

                other_cell_verts = all_cells_verts[ci]

                closest_point_dist_squared, closest_point_coords, closest_node_indices, projection_factor = calculate_closest_point_dist_squared(
                    num_vertices, this_nc, other_cell_verts, closest_ni
                )

                if (
                        closest_point_dist_squared != -1
                        and closest_point_dist_squared
                        < close_0_at
                ):
                    close_points_exist[ni][ci] = 1
                    close_points[ni][ci] = closest_point_coords
                    close_points_node_indices[ni][ci] = closest_node_indices
                    close_points_node_projection_factors[ni][
                        ci] = projection_factor
                    close_point_smoothness_factors[ni][
                        ci
                    ] = closeness_smoothening_linear_function(
                        close_0,
                        close_1,
                        np.sqrt(closest_point_dist_squared),
                    )
                elif closest_node_dist < close_0_at:
                    close_points_exist[ni][ci] = 1
                    close_points[ni][ci] = other_cell_verts[closest_ni]
                    closest_node_indices[1] = closest_ni
                    close_points_node_indices[ni][ci] = closest_node_indices
                    close_point_smoothness_factors[ni][
                        ci
                    ] = closeness_smoothening_linear_function(
                        close_0,
                        close_1,
                        np.sqrt(closest_node_dist),
                    )
                elif are_nodes_inside_other_cells[ni][ci] == 1:
                    close_points_exist[ni][ci] = 2
                    close_point_smoothness_factors[ni][ci] = 1.0
                    continue

    return (
        close_points_exist,
        close_points,
        close_points_node_indices,
        close_points_node_projection_factors,
        close_point_smoothness_factors,
    )


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def check_if_line_segment_intersects_polygon(
        a, b, normal_to_line_segment, polygon, ignore_vertex_ix
):
    num_vertices = polygon.shape[0]

    nls_x, nls_y = normal_to_line_segment

    for vi in range(num_vertices):
        next_ix = (vi + 1) % num_vertices

        if ignore_vertex_ix != -1:
            if vi == ignore_vertex_ix:
                continue
            if next_ix == ignore_vertex_ix:
                continue

        this_coords = polygon[vi]
        next_coords = polygon[next_ix]

        is_left_a = is_left(this_coords, next_coords, a)
        is_left_b = is_left(this_coords, next_coords, b)

        if is_left_a < 0 and is_left_b < 0:
            continue
        elif is_left_a > 0 and is_left_b > 0:
            continue

        if is_left_a == 0 or is_left_b == 0:
            return True

        alpha_x, alpha_y = a - this_coords
        beta_x, beta_y = next_coords - this_coords

        denominator = beta_x * nls_x + beta_y * nls_y
        if denominator == 0:
            return False

        t = (alpha_x * nls_x + alpha_y * nls_y) / denominator

        if 0 <= t < 1:
            return True
        else:
            continue

    return False


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def cross_prod(a, b):
    return (a[0] * b[1]) - (a[1] * b[0])


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def is_given_vector_between_others(x, alpha, beta):
    cp1 = cross_prod(alpha, x)
    cp2 = cross_prod(x, beta)

    if abs(cp1) < 1e-4 or abs(cp2) < 1e-4:
        return 1
    elif cp1 > 0 and cp2 > 0:
        return 1
    else:
        return 0


# -----------------------------------------------------------------


@nb.jit(nopython=True)
def check_if_line_segment_from_node_self_intersects(
        start_coord, end_coord, polygon_coords, start_coord_ix
):
    num_vertices = polygon_coords.shape[0]
    si_minus1 = (start_coord_ix - 1) % num_vertices
    si_plus1 = (start_coord_ix + 1) % num_vertices

    si_plus1_coord = polygon_coords[si_plus1]
    si_minus1_coord = polygon_coords[si_minus1]
    edge_vector_to_plus = si_plus1_coord - start_coord

    edge_vector_from_minus = start_coord - si_minus1_coord

    rough_tangent_vector = edge_vector_to_plus + edge_vector_from_minus

    ipv = rotate_vector_ccw(rough_tangent_vector)

    v = end_coord - start_coord

    if is_given_vector_between_others(v, ipv, -1 * edge_vector_from_minus) == 1:
        return 1
    if is_given_vector_between_others(v, edge_vector_to_plus, ipv) == 1:
        return 1
    else:
        return 0


# -----------------------------------------------------------------


@nb.jit(nopython=True)
def close_to_zero(x, tol):
    if x < 0:
        x *= -1

    if x < tol:
        return True

    return False


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def check_if_line_segment_intersects_vertical_line(start, end, min_y, max_y, x):
    sx, sy = start
    ex, ey = end

    if sx < x and ex < x:
        return 0
    if sx > x and ex > x:
        return 0
    else:
        # y = m*x + b
        # sy = m*sx + b
        # ey = m*ex + b
        # ey = m*ex + (sy - m*sx)
        # ey - sy = m*(ex - sx)
        # m = (ey - sy)/(ex - sx)
        # b = sy - m*sx

        denom = ex - sx
        if abs(denom) < 1e-4:
            average_x = (sx + ex) / 2.0
            if (average_x - x) < 1e-4:
                return 1
            else:
                return 0

        m = (ey - sy) / denom
        b = sy - m * sx
        y_intersect = m * x + b

        if y_intersect > max_y:
            return 0
        elif y_intersect < min_y:
            return 0
        else:
            return 1


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def check_if_line_segment_intersects_horizontal_line(start, end, min_x, max_x,
                                                     y):
    sx, sy = start
    ex, ey = end

    if sy < y and ey < y:
        return 0
    if sy > y and ey > y:
        return 0
    else:
        # y = m*x + b
        # sx = m*sy + b
        # ex = m*ey + b
        # ex = m*ey + (sx - m*sy)
        # ex - sx = m*(ey - sy)
        # m = (ex - sx)/(ey - sy)
        # b = sy - m*sx

        denom = ey - sy
        if abs(denom) < 1e-4:
            average_y = (sy + ey) / 2.0
            if (average_y - y) < 1e-4:
                return 1
            else:
                return 0

        m = (ex - sx) / denom

        b = sx - m * sy
        x_intersect = m * y + b

        if x_intersect > max_x:
            return 0
        elif x_intersect < min_x:
            return 0
        else:
            return 1


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def check_if_line_segment_intersects_box(start, end, min_x, max_x, min_y,
                                         max_y):
    if check_if_line_segment_intersects_vertical_line(start, end, min_y, max_y,
                                                      min_x):
        return 1
    if check_if_line_segment_intersects_vertical_line(start, end, min_y, max_y,
                                                      max_x):
        return 1
    if check_if_line_segment_intersects_horizontal_line(
            start, end, min_x, max_x, min_y
    ):
        return 1
    if check_if_line_segment_intersects_horizontal_line(
            start, end, min_x, max_x, max_y
    ):
        return 1

    return 0


# -----------------------------------------------------------------
# @nb.jit(nopython=True)
def check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
        pi_a,
        vi_a,
        pi_b,
        vi_b,
        all_polygon_coords,
        all_polygons_bb_coords,
):
    coords_a = all_polygon_coords[pi_a, vi_a]
    coords_b = all_polygon_coords[pi_b, vi_b]

    normal_to_line_segment = rotate_vector_ccw(coords_b - coords_a)
    check0 = check_if_line_segment_from_node_self_intersects(
        coords_a, coords_b, all_polygon_coords[pi_a], vi_a
    )
    check1 = check_if_line_segment_from_node_self_intersects(
        coords_b, coords_a, all_polygon_coords[pi_b], vi_b
    )

    if check0 == 1:
        return 100000
    elif check1 == 1:
        return 100000

    num_intersections = 0
    num_polygons = all_polygon_coords.shape[0]
    for pi in range(num_polygons):
        if pi == pi_a or pi == pi_b:
            continue
        else:
            min_x, max_x, min_y, max_y = all_polygons_bb_coords[pi]
            if (
                    check_if_line_segment_intersects_box(
                        coords_a, coords_b, min_x, max_x, min_y, max_y
                    )
                    == 1
            ):
                if (
                        check_if_line_segment_intersects_polygon(
                            coords_a,
                            coords_b,
                            normal_to_line_segment,
                            all_polygon_coords[pi],
                            -1,
                        )
                        == 1
                ):
                    num_intersections += 1

    return num_intersections


# -----------------------------------------------------------


@nb.jit(nopython=True)
def check_if_line_segment_going_from_vertex_of_one_polygon_to_point_passes_through_any_polygon(
        pi_a,
        vi_a,
        point,
        all_polygon_coords,
        all_polygons_bb_coords,
):
    coords_a = all_polygon_coords[pi_a, vi_a]

    normal_to_line_segment = rotate_vector_ccw(point - coords_a)

    #    if check_if_line_segment_from_node_self_intersects(coords_a, point, all_polygon_coords[pi_a], vi_a) == 1:
    #        return 100000
    #    else:

    num_intersections = 0
    num_polygons = all_polygon_coords.shape[0]
    for pi in range(num_polygons):
        if pi == pi_a:
            continue
        else:
            min_x, max_x, min_y, max_y = all_polygons_bb_coords[pi]
            if (
                    check_if_line_segment_intersects_box(
                        coords_a, point, min_x, max_x, min_y, max_y
                    )
                    == 1
            ):
                if (
                        check_if_line_segment_intersects_polygon(
                            coords_a,
                            point,
                            normal_to_line_segment,
                            all_polygon_coords[pi],
                            -1,
                        )
                        == 1
                ):
                    num_intersections += 1

    return num_intersections


# -----------------------------------------------------------------


@nb.jit(nopython=True)
def create_initial_line_segment_intersection_matrix(
        num_cells,
        num_vertices_per_cell,
        init_cells_bb_array,
        init_all_cells_verts,
):
    line_segment_intersection_matrix = -1 * np.ones(
        (num_cells, num_vertices_per_cell, num_cells, num_vertices_per_cell),
        dtype=np.int64,
    )
    info_update_tracker = np.zeros_like(
        line_segment_intersection_matrix, dtype=np.int64
    )

    for ci in range(num_cells):
        for ni in range(num_vertices_per_cell):
            for other_ci in range(num_cells):
                if ci != other_ci:
                    relevant_info_update_tracker_slice = \
                        info_update_tracker[ci][ni][
                            other_ci
                        ]
                    for other_ni in range(num_vertices_per_cell):
                        if relevant_info_update_tracker_slice[other_ni] != 1:
                            info_update_tracker[ci][ni][other_ci][other_ni] = 1
                            info_update_tracker[other_ci][other_ni][ci][ni] = 1

                            does_line_segment_between_nodes_intersect = check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
                                ci,
                                ni,
                                other_ci,
                                other_ni,
                                init_all_cells_verts,
                                init_cells_bb_array,
                            )

                            line_segment_intersection_matrix[ci][ni][other_ci][
                                other_ni
                            ] = does_line_segment_between_nodes_intersect
                            line_segment_intersection_matrix[other_ci][
                                other_ni][ci][
                                ni
                            ] = does_line_segment_between_nodes_intersect

    return line_segment_intersection_matrix


# -----------------------------------------------------------------


@nb.jit(nopython=True)
def update_line_segment_intersection_matrix(
        last_updated_cell_ix,
        num_cells,
        num_vertices_per_cell,
        all_cells_verts,
        cells_bb_array,
        line_segment_intersection_matrix,
):
    # 1 if info has been updated, 0 otherwise
    info_update_tracker = np.zeros_like(
        line_segment_intersection_matrix, dtype=np.int64
    )

    for ni in range(num_vertices_per_cell):
        for other_ci in range(num_cells):
            if other_ci != last_updated_cell_ix:
                relevant_info_update_tracker_slice = info_update_tracker[
                    last_updated_cell_ix
                ][ni][other_ci]
                for other_ni in range(num_vertices_per_cell):
                    if relevant_info_update_tracker_slice[other_ni] != 1:
                        info_update_tracker[last_updated_cell_ix][ni][other_ci][
                            other_ni
                        ] = 1
                        info_update_tracker[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = 1

                        does_line_segment_between_nodes_intersect = check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
                            last_updated_cell_ix,
                            ni,
                            other_ci,
                            other_ni,
                            all_cells_verts,
                            cells_bb_array,
                        )

                        line_segment_intersection_matrix[last_updated_cell_ix][
                            ni][
                            other_ci
                        ][other_ni] = does_line_segment_between_nodes_intersect
                        line_segment_intersection_matrix[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = does_line_segment_between_nodes_intersect

    return line_segment_intersection_matrix


# -----------------------------------------------------------------


# @nb.jit(nopython=True)
def init_lseg_intersects_and_dist_sq_matrices_old(
        num_cells,
        num_vertices_per_cell,
        init_cells_bb_array,
        init_all_cells_verts,
):
    distance_squared_matrix = -1 * np.ones(
        (num_cells, num_vertices_per_cell, num_cells, num_vertices_per_cell),
        dtype=np.float64,
    )
    line_segment_intersection_matrix = -1 * np.ones(
        (num_cells, num_vertices_per_cell, num_cells, num_vertices_per_cell),
        dtype=np.int64,
    )

    info_update_tracker = np.zeros_like(distance_squared_matrix, dtype=np.int64)

    for ci in range(num_cells):
        for ni in range(num_vertices_per_cell):
            for other_ci in range(num_cells):
                if ci != other_ci:
                    relevant_info_update_tracker_slice = \
                        info_update_tracker[ci][ni][
                            other_ci
                        ]
                    for other_ni in range(num_vertices_per_cell):
                        if relevant_info_update_tracker_slice[other_ni] != 1:
                            info_update_tracker[ci][ni][other_ci][other_ni] = 1
                            info_update_tracker[other_ci][other_ni][ci][ni] = 1

                            does_line_segment_between_nodes_intersect = check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
                                ci,
                                ni,
                                other_ci,
                                other_ni,
                                init_all_cells_verts,
                                init_cells_bb_array,
                            )

                            line_segment_intersection_matrix[ci][ni][other_ci][
                                other_ni
                            ] = does_line_segment_between_nodes_intersect
                            line_segment_intersection_matrix[other_ci][
                                other_ni][ci][
                                ni
                            ] = does_line_segment_between_nodes_intersect

                            this_node = init_all_cells_verts[ci][ni]
                            other_node = init_all_cells_verts[other_ci][
                                other_ni]

                            squared_dist = calculate_squared_dist(this_node,
                                                                  other_node)
                            distance_squared_matrix[ci][ni][other_ci][
                                other_ni
                            ] = squared_dist
                            distance_squared_matrix[other_ci][other_ni][ci][
                                ni
                            ] = squared_dist

    return distance_squared_matrix, line_segment_intersection_matrix


# -----------------------------------------------------------------


# @nb.jit(nopython=True, nogil=True)
def dist_squared_and_line_segment_calculation_worker(
        dist_squared_matrix,
        line_segment_intersect_matrix,
        polygon_bbs,
        polygons,
        task_addresses,
):
    for n in range(task_addresses.shape[0]):
        a_ci, a_ni, b_ci, b_ni = (
            task_addresses[n][0],
            task_addresses[n][1],
            task_addresses[n][2],
            task_addresses[n][3],
        )
        dist_squared = calculate_squared_dist(
            polygons[a_ci][a_ni], polygons[b_ci][b_ni]
        )
        num_intersects = check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
            a_ci,
            a_ni,
            b_ci,
            b_ni,
            polygons,
            polygon_bbs,
        )

        dist_squared_matrix[a_ci][a_ni][b_ci][b_ni] = dist_squared
        dist_squared_matrix[b_ci][b_ni][a_ci][a_ni] = dist_squared
        line_segment_intersect_matrix[a_ci][a_ni][b_ci][b_ni] = num_intersects
        line_segment_intersect_matrix[b_ci][b_ni][a_ci][a_ni] = num_intersects


@nb.jit("void(float64[:,:,:,:], float64[:,:,:], int64[:,:])", nopython=True,
        nogil=True)
def dist_squared_calculation_worker(dist_squared_matrix, polygons,
                                    task_addresses):
    for n in range(task_addresses.shape[0]):
        a_ci, a_ni, b_ci, b_ni = (
            task_addresses[n][0],
            task_addresses[n][1],
            task_addresses[n][2],
            task_addresses[n][3],
        )
        dist_squared = calculate_squared_dist(
            polygons[a_ci][a_ni], polygons[b_ci][b_ni]
        )

        dist_squared_matrix[a_ci][a_ni][b_ci][b_ni] = dist_squared
        dist_squared_matrix[b_ci][b_ni][a_ci][a_ni] = dist_squared


@nb.jit(nopython=True)
def create_dist_and_line_segment_intersection_test_args_relative_to_specific_cell(
        specific_cell_ix, num_cells, num_vertices_per_cell
):
    tasks = []

    for ni in range(num_vertices_per_cell):
        for other_ci in range(num_cells):
            if other_ci != specific_cell_ix:
                for other_ni in range(num_vertices_per_cell):
                    tasks.append((specific_cell_ix, ni, other_ci, other_ni))

    return tasks


@nb.jit(nopython=True, nogil=True)
def create_dist_and_line_segment_intersection_test_args(
        num_cells, num_vertices_per_cell
):
    info_update_tracker = -1 * np.ones(
        (num_cells, num_vertices_per_cell, num_cells, num_vertices_per_cell),
        dtype=np.float64,
    )
    tasks = []

    for ci in range(num_cells):
        for ni in range(num_vertices_per_cell):
            for other_ci in range(num_cells):
                if ci != other_ci:
                    relevant_info_update_tracker_slice = \
                        info_update_tracker[ci][ni][
                            other_ci
                        ]
                    for other_ni in range(num_vertices_per_cell):
                        if relevant_info_update_tracker_slice[other_ni] != 1:
                            info_update_tracker[ci][ni][other_ci][other_ni] = 1
                            info_update_tracker[other_ci][other_ni][ci][ni] = 1
                            tasks.append((ci, ni, other_ci, other_ni))

    return tasks


def create_initial_line_segment_intersection_and_dist_squared_matrices(
        num_threads,
        tasks,
        num_cells,
        num_vertices,
        polygon_bbs,
        polygons,
        sequential=True,
):
    num_tasks = tasks.shape[0]

    dist_squared_matrix = -1 * np.ones(
        (num_cells, num_vertices, num_cells, num_vertices), dtype=np.float64
    )
    line_segment_intersect_matrix = -1 * np.ones(
        (num_cells, num_vertices, num_cells, num_vertices), dtype=np.int64
    )

    if num_tasks != 0:
        chunklen = (num_tasks + num_threads - 1) // num_threads
        # Create argument tuples for each input chunk
        chunks = []
        for i in range(num_threads):
            relevant_tasks = tasks[i * chunklen: (i + 1) * chunklen]
            chunks.append(
                (
                    dist_squared_matrix,
                    line_segment_intersect_matrix,
                    polygon_bbs,
                    polygons,
                    relevant_tasks,
                )
            )

        # Spawn one thread per chunk
        if not sequential:
            threads = [
                threading.Thread(
                    target=dist_squared_and_line_segment_calculation_worker,
                    args=c
                )
                for c in chunks
            ]

            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        else:
            for chunk in chunks:
                dist_squared_and_line_segment_calculation_worker(*chunk)

    return dist_squared_matrix, line_segment_intersect_matrix


# -----------------------------------------------------------------
@nb.jit(nopython=True)
def update_dist_squared_matrix(
        last_updated_cell_ix,
        num_cells,
        num_vertices_per_cell,
        all_cells_verts,
        distance_squared_matrix,
):
    # 1 if info has been updated, 0 otherwise
    info_update_tracker = np.zeros_like(distance_squared_matrix, dtype=np.int64)

    for ni in range(num_vertices_per_cell):
        for other_ci in range(num_cells):
            if other_ci != last_updated_cell_ix:
                relevant_info_update_tracker_slice = info_update_tracker[
                    last_updated_cell_ix
                ][ni][other_ci]
                for other_ni in range(num_vertices_per_cell):
                    if relevant_info_update_tracker_slice[other_ni] != 1:
                        info_update_tracker[last_updated_cell_ix][ni][other_ci][
                            other_ni
                        ] = 1
                        info_update_tracker[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = 1

                        this_node = all_cells_verts[last_updated_cell_ix][ni]
                        other_node = all_cells_verts[other_ci][other_ni]
                        squared_dist = calculate_squared_dist(this_node,
                                                              other_node)

                        distance_squared_matrix[last_updated_cell_ix][ni][
                            other_ci][
                            other_ni
                        ] = squared_dist
                        distance_squared_matrix[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = squared_dist

    return distance_squared_matrix


# ---------------------------------------------------------------------------------
@nb.jit(nopython=True)
def update_line_segment_intersection_and_dist_squared_matrices_old(
        last_updated_cell_ix,
        num_cells,
        num_vertices_per_cell,
        all_cells_verts,
        cells_bb_array,
        distance_squared_matrix,
        line_segment_intersection_matrix,
):
    # 1 if info has been updated, 0 otherwise
    info_update_tracker = np.zeros_like(
        line_segment_intersection_matrix, dtype=np.int64
    )

    for ni in range(num_vertices_per_cell):
        for other_ci in range(num_cells):
            pass
            if other_ci != last_updated_cell_ix:
                relevant_info_update_tracker_slice = info_update_tracker[
                    last_updated_cell_ix
                ][ni][other_ci]
                for other_ni in range(num_vertices_per_cell):
                    if relevant_info_update_tracker_slice[other_ni] != 1:
                        info_update_tracker[last_updated_cell_ix][ni][other_ci][
                            other_ni
                        ] = 1
                        info_update_tracker[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = 1

                        does_line_segment_between_nodes_intersect = check_if_line_segment_going_from_vertex_of_one_polygon_to_vertex_of_another_passes_through_any_polygon(
                            last_updated_cell_ix,
                            ni,
                            other_ci,
                            other_ni,
                            all_cells_verts,
                            cells_bb_array,
                        )

                        line_segment_intersection_matrix[last_updated_cell_ix][
                            ni][
                            other_ci
                        ][other_ni] = does_line_segment_between_nodes_intersect
                        line_segment_intersection_matrix[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = does_line_segment_between_nodes_intersect

                        this_node = all_cells_verts[last_updated_cell_ix][ni]
                        other_node = all_cells_verts[other_ci][other_ni]
                        squared_dist = calculate_squared_dist(this_node,
                                                              other_node)

                        distance_squared_matrix[last_updated_cell_ix][ni][
                            other_ci][
                            other_ni
                        ] = squared_dist
                        distance_squared_matrix[other_ci][other_ni][
                            last_updated_cell_ix
                        ][ni] = squared_dist

    return distance_squared_matrix, line_segment_intersection_matrix


def update_line_segment_intersection_and_dist_squared_matrices(
        num_threads,
        given_tasks,
        all_cells_verts,
        cells_bb_array,
        distance_squared_matrix,
        line_segment_intersection_matrix,
        sequential=False,
):
    num_tasks = given_tasks.shape[0]

    if num_tasks != 0:
        chunklen = (num_tasks + num_threads - 1) // num_threads
        # Create argument tuples for each input chunk
        chunks = []
        for i in range(num_threads):
            relevant_tasks = given_tasks[i * chunklen: (i + 1) * chunklen]
            chunks.append(
                (
                    distance_squared_matrix,
                    line_segment_intersection_matrix,
                    cells_bb_array,
                    all_cells_verts,
                    relevant_tasks,
                )
            )

        # Spawn one thread per chunk
        if not sequential:
            threads = [
                threading.Thread(
                    target=dist_squared_and_line_segment_calculation_worker,
                    args=c
                )
                for c in chunks
            ]

            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        else:
            for chunk in chunks:
                dist_squared_and_line_segment_calculation_worker(*chunk)

    return distance_squared_matrix, line_segment_intersection_matrix


if __name__ == "__main__":
    print("=============== geometry.py ====================")
    polygon = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=np.float64)
    unit_inside_pointing_vecs = calculate_unit_inside_pointing_vecs(polygon)
    print("================================================")
