import copy
import py_model.hardio as hio
import numpy as np


def is_pos(label):
    return label == "poly"


def is_force(label):
    return np.any([(x == label) for x in ["rgtp_forces", "cyto_forces",
                                          "sum_forces", "edge_forces",
                                          "edge_forces_minus"]])


def is_other_vector(label):
    return np.any([(x in label) for x in ["uivs, edge_vecs, uevs, uevs_minus"]])


def is_vector(label):
    return is_force(label) or is_other_vector(label)


def is_extra(label):
    return is_pos_extra(label) or is_force(label) or \
           is_other_vector(label)


def flatten(items):
    flattened = []
    for x in items:
        if type(x) == list or type(x) == tuple:
            for y in x:
                flattened.append(y)
        else:
            flattened.append(x)
    return flattened


LABEL_COLOR_DICT = dict()
# BASICS = ["poly", "rac_acts", "rac_inacts", "rho_acts", "rho_inacts",
#           "sum_forces"]
BASIC_COLORS = ["black", "blue", "seagreen", "red", "darkviolet", "orangered"]
LABEL_COLOR_DICT.update(zip(hio.BASICS, BASIC_COLORS))
# GEOMETRY = ["uivs"]
GEOMETRY_COLORS = ["mediumseagreen"]
LABEL_COLOR_DICT.update(zip(hio.GEOMETRY, GEOMETRY_COLORS))
# RAC_RATES = ["kgtps_rac", "kdgtps_rac"]
RAC_RATE_COLORS = ["blue", "seagreen"]
LABEL_COLOR_DICT.update(zip(hio.RAC_RATES, RAC_RATE_COLORS))
# RHO_RATES = ["kgtps_rho", "kdgtps_rho"]
RHO_RATE_COLORS = ["red", "darkviolet"]
LABEL_COLOR_DICT.update(zip(hio.RHO_RATES, RHO_RATE_COLORS))
# FORCES = ["rgtp_forces", "edge_forces", "cyto_forces"]
FORCE_COLORS = ["darkslategrey", "purple", "darkgoldenrod"]
LABEL_COLOR_DICT.update(zip(hio.FORCES, FORCE_COLORS))
# CALC_KGTPS_RAC = ["x_cils", "x_coas"]
CALC_KGTPS_RAC_COLORS = ["indianred", "forestgreen"]
LABEL_COLOR_DICT.update(zip(hio.CALC_KGTPS_RAC, CALC_KGTPS_RAC_COLORS))
# OTHERS = ["edge_strains", "uevs", "poly_area", "coa_update", "cil_update"]
OTHER_COLORS = ["purple", "skyblue", "tomato", "forestgreen", "indianred"]
LABEL_COLOR_DICT.update(zip(hio.OTHERS, OTHER_COLORS))

EXTRA_POSITIONS = \
    [("poly_x", "poly_y"), "centroid"]
EXTRA_POSITION_COLORS = [("dimgrey", "darkgrey"), "rosybrown"]
LABEL_COLOR_DICT.update(
    zip(flatten(EXTRA_POSITIONS), flatten(EXTRA_POSITION_COLORS)))

EXTRA_FORCES = [
    ("{}_x".format(label), "{}_y".format(label))
    for label in ["sum_forces"] + hio.FORCES
]
EXTRA_FORCE_COLORS = [("tomato", "lightcoral")] + \
                     [("cadetblue", "powderblue"), ("violet", "thistle")] + \
                     [("yellowgreen", "palegreen"), ("gold", "wheat")]
LABEL_COLOR_DICT.update(zip(flatten(EXTRA_FORCES), flatten(EXTRA_FORCE_COLORS)))

EXTRA_LABELS = flatten(EXTRA_POSITIONS) + flatten(EXTRA_FORCES)
ALL_LABELS = hio.BASICS + hio.FORCES + EXTRA_POSITIONS + EXTRA_FORCES + \
             hio.GEOMETRY + hio.RAC_RATES + hio.RHO_RATES + \
             hio.CALC_KGTPS_RAC + hio.OTHERS
ALL_LABELS_NO_EXTRAS = hio.BASICS + hio.FORCES + \
                       hio.GEOMETRY + hio.RAC_RATES + hio.RHO_RATES + \
                       hio.CALC_KGTPS_RAC + hio.OTHERS
VECTOR_LABELS = ["poly"] + ["sum_forces"] + hio.FORCES


def calc_min_maxes_for_labels(labels, py_data, rust_data, max_step):
    num_cells = len(py_data)

    all_dat = np.zeros(0, dtype=np.float64)
    for label in labels:
        for ix in range(num_cells):
            pycd = py_data[ix]
            rucd = rust_data[ix]
            all_dat = np.append(all_dat, pycd[label][:max_step].flatten())
            all_dat = np.append(all_dat, rucd[label][:max_step].flatten())

    min_lim = np.min([0.0, 1.2 * np.min(all_dat)])
    max_lim = 1.2 * np.max(all_dat)
    if abs(min_lim - max_lim) < 1e-8:
        min_lim = max_lim - 0.5
        max_lim += 0.5
    return min_lim, max_lim


def calc_min_maxes_given_label_grouping(label_groups, py_data, rust_data,
                                        max_step):
    min_maxes = []
    for label_group in label_groups:
        if type(label_group) != tuple:
            label = label_group
            min_maxes.append((label, calc_min_maxes_for_labels([label],
                                                               py_data,
                                                               rust_data,
                                                               max_step)))
        elif type(label_group) == tuple:
            group_min_max = calc_min_maxes_for_labels(label_group, py_data,
                                                      rust_data, max_step)
            min_maxes.append((label_group, group_min_max))
        else:
            raise Exception("received a label grouping that is not packed in "
                            "a tuple: {}".format(label_group))
    return dict(min_maxes)


class PlotDataGroup:
    labels = []
    label_groups = []
    description = ""

    def __init__(self, py_data_dict, rust_data_dict, labels, label_groups,
                 description):
        self.labels = labels
        self.label_groups = label_groups
        self.description = description
        self.py_dat = [dict([(label, py_data_dict[ix][label])
                             for label in self.labels]) for ix
                       in range(len(py_data_dict))]
        self.rust_dat = [dict([(label, rust_data_dict[ix][label])
                               for label in self.labels]) for ix
                         in range(len(rust_data_dict))]
        self.grouped_ylims_dict = dict()
        self.ylims_dict = dict()
        self.recalc_ylims(np.min([len(self.py_dat[0]), len(self.rust_dat[0])]))

    def recalc_ylims(self, max_step):
        self.grouped_ylims_dict = calc_min_maxes_given_label_grouping(
            self.label_groups, self.py_dat, self.rust_dat, max_step)
        self.ylims_dict = calc_min_maxes_given_label_grouping(
            self.labels, self.py_dat, self.rust_dat, max_step)


def split_vector_dat_by_components(label, dat_per_step):
    if "_x" in label:
        return copy.deepcopy(dat_per_step[:, :, 0])
    elif "_y" in label:
        return copy.deepcopy(dat_per_step[:, :, 1])
    else:
        raise Exception("Could not extract as component: {}.".format(label))


def get_comp_ix(label):
    if label[-2:] == "_x":
        return label[:-2], 0
    else:
        return label[:-2], 1


def calc_extras(extra_labels, data_dict):
    cells = np.arange(len(data_dict))
    for ci in cells:
        data_dict[ci]["uevs_minus"] = np.roll(data_dict[ci]["uevs"], -1, axis=1)
    for extra in extra_labels:
        for ci in cells:
            if extra == "centroid":
                poly_per_step = data_dict[ci]["poly"]
                centroids = np.average(poly_per_step, axis=1)
                centroid_dists = np.linalg.norm(centroids, axis=1)
                data_dict[ci][extra] = centroid_dists
            elif extra[-2:] in ["_x", "_y"]:
                parent, comp_ix = get_comp_ix(extra)
                data_dict[ci][extra] = data_dict[ci][parent][:, :, comp_ix]
    return data_dict


def norm_vector_data(labels, data_dict):
    for label in labels:
        for ci in range(len(data_dict)):
            data_dict[ci][label] = np.linalg.norm(data_dict[ci][label], axis=2)
    return data_dict


def prep_data(py_data_dict_per_cell, rust_data_dict_per_cell):
    py_data_dict_per_cell = calc_extras(EXTRA_LABELS, py_data_dict_per_cell)
    rust_data_dict_per_cell = calc_extras(EXTRA_LABELS, rust_data_dict_per_cell)
    py_data_dict_per_cell = norm_vector_data(VECTOR_LABELS,
                                             py_data_dict_per_cell)
    rust_data_dict_per_cell = norm_vector_data(VECTOR_LABELS,
                                               rust_data_dict_per_cell)

    return py_data_dict_per_cell, rust_data_dict_per_cell
