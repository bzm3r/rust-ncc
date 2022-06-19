import os
import copy
import orjson
import numpy as np

HEADER_LABELS = ["num_cells", "num_int_steps", "vertex_eta",
                 "f", "l", "t", "l3d", "kgtp", "zero_at",
                 "one_at", "cil_mag", "coa_los_penalty",
                 "coa_halfmax_dist",
                 "coa_vertex_mag", "coa_distrib_exp", "vertex_eta", "cell_r",
                 "rest_edge_len", "rest_area", "stiffness_edge",
                 "const_protrusive", "const_retractive", "stiffness_cyto",
                 "k_mem_on_vertex", "k_mem_off", "diffusion_rgtp", "init_rac",
                 "init_rho", "halfmax_vertex_rgtp",
                 "halfmax_vertex_rgtp_conc",
                 "kgtp_rac", "kgtp_rac_auto", "kdgtp_rac", "kdgtp_rho_on_rac",
                 "halfmax_tension_inhib", "tension_inhib", "kgtp_rho",
                 "kgtp_rho_auto", "kdgtp_rho", "kdgtp_rac_on_rho",
                 "randomization", "rand_avg_t", "rand_std_t", "rand_mag",
                 "num_rand_vs"]
BASICS = ["tpoint", "poly", "rac_acts", "rac_inacts", "rho_acts",
          "rho_inacts", "sum_forces"]
GEOMETRY = ["uivs"]
FORCES = ["rgtp_forces", "edge_forces", "cyto_forces"]
RAC_RATES = ["kgtps_rac", "kdgtps_rac"]
RHO_RATES = ["kgtps_rho", "kdgtps_rho"]
CALC_KGTPS_RAC = ["x_cils", "x_coas"]
DIFFUSION = ["rac_act_net_fluxes", "rac_inact_net_fluxes",
             "rho_act_net_fluxes", "rho_inact_net_fluxes"]
OTHERS = ["edge_strains", "avg_tens_strain", "poly_area"] + DIFFUSION + [
    "x_tens"]

DATA_LABELS = BASICS + GEOMETRY + FORCES + RAC_RATES + RHO_RATES + \
              CALC_KGTPS_RAC + OTHERS
FILE_NAME_TEMPLATE = "out_euler_T={final_t}_E={num_int_steps}_NC={" \
                     "num_cells}_CIL={cil}_COA={coa}.{ext}"

NO_LOG = 100
BASIC_INFO = 99
DYNAMICS_INFO = 98
COA_INFO = 97
VOL_EX_INFO = 96
DETAILED_VOL_EX_INFO = 95
GEOM_ALG_INFO = 94


def make_params_list(params):
    params_list = []
    for req_label in HEADER_LABELS:
        data = params[req_label]
        if type(data) == np.float64:
            store = float(data)
        elif type(data) == np.ndarray:
            store = [float(x) for x in data]
        else:
            store = data
        params_list.append((req_label, store))
    return params_list


def validate_data(req_labels, data):
    for req_label, d in zip(req_labels, data):
        label, value = d
        if req_label != label:
            raise Exception("Expected label {}, found {}".format(req_label,
                                                                 label))


def validate_buffer(expected_length, buf, buf_description):
    if len(buf) != expected_length:
        raise Exception("{} expected {} elements. Found: {}".format(
            buf_description, expected_length, len(buf)))


class Writer:
    def __init__(self, out_dir, name, log_level, params, snap_period):
        global HEADER_LABELS
        global BASICS
        global GEOMETRY
        global RAC_RATES
        global RHO_RATES
        global FORCES
        global CALC_KGTPS_RAC
        global OTHERS
        self.header_labels = copy.deepcopy(HEADER_LABELS)
        self.basics = copy.deepcopy(BASICS)
        self.geometry = copy.deepcopy(GEOMETRY)
        self.rac_rates = copy.deepcopy(RAC_RATES)
        self.rho_rates = copy.deepcopy(RHO_RATES)
        self.forces = copy.deepcopy(FORCES)
        self.calc_kgtps_rac = copy.deepcopy(CALC_KGTPS_RAC)
        self.others = copy.deepcopy(OTHERS)
        self.data_labels = self.basics + self.geometry + self.forces + \
                           self.rac_rates + self.rho_rates + \
                           self.calc_kgtps_rac + self.others

        self.params = params
        self.num_int_steps = self.params["num_int_steps"]
        self.num_cells = self.params["num_cells"]
        self.cil_mag = self.params["cil_mag"]
        self.coa_mag = self.params["coa_mag"]

        self.out_dir = out_dir
        self.name = name
        self.write_file_path = os.path.join(self.out_dir, self.name + ".dat")
        self.log_file_path = os.path.join(self.out_dir, self.name + ".log")
        print("saving to: {}".format(self.write_file_path))
        # logging.basicConfig(filename=self.log_file_path, filemode="w",
        #                     format="%(message)s", level=log_level)
        self.curr_int_step_ix = 0
        self.can_save_int_step = False
        self.can_save_initial = False
        self.curr_cell_data = {"int_steps": []}
        self.curr_cell_ix = 0
        self.can_save_cell = False
        self.curr_step_data = {"cells": []}
        self.data = dict([("header", {}), ("tpoints", [])])
        self.finished = False
        self.last_saved = None
        self.snap_period = snap_period
        self.snap_save = False

    def begin_cell_save(self):
        if self.snap_save:
            if self.curr_cell_ix < self.num_cells and self.curr_int_step_ix == 0:
                self.can_save_int_step = True
            else:
                msg = "cannot begin cell save: " + \
                      "curr_cell = {} (max: {}), " + \
                      "curr_int_step = {} (max: {})"
                raise Exception(msg.format(self.curr_cell_ix, self.num_cells,
                                           self.curr_int_step, self.num_int_steps))

    def begin_cell_initial_save(self):
        if self.curr_cell_ix < self.num_cells and self.curr_int_step_ix == 0:
            self.can_save_initial = True
        else:
            msg = "cannot begin cell save: " + \
                  "curr_cell = {} (max: {}), " + \
                  "curr_int_step = {} (max: {})"
            raise Exception(msg.format(self.curr_cell_ix, self.num_cells,
                                       self.curr_int_step, self.num_int_steps))

    def save_int_step(self, data):
        if self.snap_save:
            if self.can_save_int_step:
                if len(self.curr_cell_data["int_steps"]) < self.num_int_steps:
                    validate_data(self.data_labels, data)
                    self.curr_cell_data["int_steps"].append(dict(data))
                    self.curr_int_step_ix += 1
                else:
                    msg = "cannot save int step, " + \
                          "already have {} int steps stored (max: {})"
                    raise Exception(msg.format(
                        len(self.curr_cell_data["int_steps"]), self.num_int_steps)
                    )
            else:
                raise Exception("can_save_int_step is false")

    def save_initial_cell(self, data):
        if self.can_save_initial and self.curr_int_step_ix == 0:
            validate_data(self.data_labels, data)
            self.curr_cell_data["int_steps"].append(dict(data))
            self.curr_int_step_ix += 1
        else:
            raise Exception("can_save_int_step is false")

    def finish_cell_save(self):
        if self.snap_save:
            if self.curr_int_step_ix == self.num_int_steps:
                self.can_save_int_step = False
                if self.can_save_cell:
                    if len(self.curr_step_data["cells"]) < self.num_cells:
                        self.curr_step_data["cells"].append(
                            copy.deepcopy(self.curr_cell_data)
                        )
                        self.curr_cell_ix += 1
                    else:
                        msg = "cannot save save cell, " + \
                              "already have {} cells stored (max: {})"
                        raise Exception(msg.format(len(self.curr_step_data["cells"]),
                                                   self.num_cells))
                else:
                    raise Exception("can_save_cell is false")
                self.curr_cell_data = {"int_steps": []}
                self.curr_int_step_ix = 0
            else:
                msg = "cannot finish cell save; only have {} int steps " + \
                      "(required: {})"
                raise Exception(msg.format(self.curr_int_step_ix, self.num_int_steps))

    def finish_cell_initial_save(self):
        if self.can_save_cell:
            if len(self.curr_step_data["cells"]) <= self.num_cells:
                self.curr_step_data["cells"].append(
                    copy.deepcopy(self.curr_cell_data)
                )
                self.curr_cell_ix += 1
            else:
                msg = "cannot save save cell, " + \
                      "already have {} cells stored (max: {})"
                raise Exception(msg.format(len(self.curr_step_data["cells"]),
                                           self.num_cells))
        else:
            raise Exception("can_save_cell is false")
        self.curr_cell_data = {"int_steps": []}
        self.curr_int_step_ix = 0

    def begin_tpoint_save(self, curr_tpoint):
        if self.last_saved is None or not (curr_tpoint - self.last_saved <
                                           self.snap_period):
            self.snap_save = True
            self.last_saved = curr_tpoint
            if self.curr_cell_ix == 0 and self.curr_int_step_ix == 0:
                self.can_save_cell = True
            else:
                msg = "cannot begin cell save: " + \
                      "curr_cell = {} (max: {})," + \
                      "curr_int_step = {} (max: {})"
                raise Exception(
                    msg.format(
                        self.curr_cell_ix, self.num_cells,
                        self.curr_int_step_ix, self.num_int_steps)
                )

    def begin_initial_save(self):
        if self.curr_cell_ix == 0 and self.curr_int_step_ix == 0:
            self.can_save_cell = True
        else:
            msg = "cannot begin cell save: curr_cell = {} (max: {}), " + \
                  "curr_int_step = {} (max: {})"
            raise Exception(msg.format(self.curr_cell_ix, self.num_cells,
                                       self.curr_int_step_ix,
                                       self.num_int_steps))

    def finish_tpoint_save(self):
        if self.snap_save:
            if self.curr_cell_ix == self.num_cells:
                if self.can_save_int_step or self.curr_int_step_ix > 0:
                    msg = "cannot finish tstep: int step saving has not been " + \
                          "finished (can still save more int steps in current cell)"
                    raise Exception(msg)
                else:
                    self.can_save_cell = False
                    self.curr_cell_ix = 0
                    self.data["tpoints"].append(copy.deepcopy(self.curr_step_data))
                    self.curr_step_data = {"cells": []}
                    self.snap_save = False
            else:
                msg = "cannot finish tstep save; only have {} cells (required: {})"
                raise Exception(msg.format(self.curr_cell_ix, self.num_cells))

    def finish_initial_save(self):
        if self.curr_cell_ix == self.num_cells:
            if self.can_save_int_step or self.curr_int_step_ix > 0:
                msg = "cannot finish tstep: int step saving has not been " + \
                      "finished (can still save more int steps in current cell)"
                raise Exception(msg)
            else:
                self.can_save_cell = False
                self.curr_cell_ix = 0
                self.data["tpoints"].append(copy.deepcopy(self.curr_step_data))
                self.curr_step_data = {"cells": []}
        else:
            msg = "cannot finish tstep save; only have {} cells (required: {})"
            raise Exception(msg.format(self.curr_cell_ix, self.num_cells))

    def save_header(self, data):
        data = make_params_list(data)
        validate_data(self.header_labels, data)
        self.data["header"] = dict(data)

    def finish(self):
        if self.can_save_int_step:
            raise Exception("cannot finish; int step saving has not been "
                            "finished (can still save more int steps into "
                            "current cell)")
        if self.can_save_cell:
            raise Exception("cannot finish; tstep saving has not been finished "
                            "(can still save more cells into current step)")

        if os.path.exists(self.write_file_path):
            os.remove(self.write_file_path)
        with open(self.write_file_path, "wb") as f:
            f.write(orjson.dumps(self.data))
        # logging.shutdown()
