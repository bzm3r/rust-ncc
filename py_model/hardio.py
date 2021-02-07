import os
import copy
import orjson
import numpy as np

HEADER_LABELS = ["num_tsteps", "num_cells", "num_int_steps", "eta",
                 "f", "l", "t", "l3d", "kgtp", "kdgtp", "close_zero_at",
                 "close_one_at", "cil_mag", "coa_los_penalty", "coa_range",
                 "coa_mag", "coa_distrib_exp", "vertex_eta", "cell_r",
                 "rest_edge_len", "rest_area", "stiffness_edge",
                 "const_protrusive", "const_retractive", "stiffness_cyto",
                 "k_mem_on_vertex", "k_mem_off", "diffusion_rgtp", "init_rac",
                 "init_rho", "halfmax_vertex_rgtp_act",
                 "halfmax_vertex_rgtp_conc", "tot_rac", "tot_rho",
                 "kgtp_rac", "kgtp_rac_auto", "kdgtp_rac", "kdgtp_rho_on_rac",
                 "halfmax_tension_inhib", "tension_inhib", "kgtp_rho",
                 "kgtp_rho_auto", "kdgtp_rho", "kdgtp_rac_on_rho",
                 "randomization", "rand_avg_t", "rand_std_t", "rand_mag",
                 "num_rand_vs", "total_rgtp"]
BASICS = ["poly", "rac_acts", "rac_inacts", "rho_acts", "rho_inacts",
          "sum_forces"]
GEOMETRY = ["uivs"]
RAC_RATES = ["kgtps_rac", "kdgtps_rac"]
RHO_RATES = ["kgtps_rho", "kdgtps_rho"]
FORCES = ["rgtp_forces", "edge_forces", "edge_forces_minus", "cyto_forces"]
CALC_KGTPS_RAC = ["conc_rac_acts", "x_cils", "x_coas"]
DIFFUSION = ["rac_act_net_fluxes"]
OTHERS = ["edge_strains", "uevs", "poly_area", "coa_update", "cil_update"]

DATA_LABELS = BASICS + GEOMETRY + RAC_RATES + RHO_RATES + FORCES + \
              CALC_KGTPS_RAC + DIFFUSION + OTHERS

WRITE_FOLDER = "B:\\rust-ncc\\model-comparison\\py-out\\"
WRITE_FILE_NAME_TEMPLATE = "out_euler_T={}_E={}_NC={}_CIL={}_COA={}.dat"


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
    def __init__(self, params):
        global HEADER_LABELS
        global BASICS
        global GEOMETRY
        global RAC_RATES
        global RHO_RATES
        global FORCES
        global CALC_KGTPS_RAC
        global DIFFUSION
        global OTHERS
        global WRITE_FOLDER
        global WRITE_FILE_NAME_TEMPLATE
        self.header_labels = copy.deepcopy(HEADER_LABELS)
        self.basics = copy.deepcopy(BASICS)
        self.geometry = copy.deepcopy(GEOMETRY)
        self.rac_rates = copy.deepcopy(RAC_RATES)
        self.rho_rates = copy.deepcopy(RHO_RATES)
        self.forces = copy.deepcopy(FORCES)
        self.calc_kgtps_rac = copy.deepcopy(CALC_KGTPS_RAC)
        self.diffusion = copy.deepcopy(DIFFUSION)
        self.others = copy.deepcopy(OTHERS)
        self.data_labels = self.basics + self.geometry + self.forces + \
                           self.rac_rates + self.rho_rates + \
                           self.calc_kgtps_rac + self.diffusion + self.others

        self.params = params

        self.write_folder = WRITE_FOLDER
        self.write_file_name_template = WRITE_FILE_NAME_TEMPLATE
        self.write_file_path_template = self.write_folder + \
                                        self.write_file_name_template
        self.write_file_path = ""
        self.write_file_label = ""

        self.num_tsteps = self.params["num_tsteps"]
        self.num_int_steps = self.params["num_int_steps"]
        self.num_cells = self.params["num_cells"]
        self.cil_mag = self.params["cil_mag"]
        self.coa_mag = self.params["coa_mag"]

        self.int_step_buffer = []
        self.current_cell_ix = None
        self.tstep_buffer = [list() for _ in range(self.num_cells)]
        self.main_buffer = dict([("header", {}), ("tsteps", [])])
        self.finished = False

        self.__init_writer()

    def __init_writer(self):
        self.write_file_path = self.write_file_path_template.format(
            self.params["num_tsteps"], self.params["num_int_steps"],
            self.params["num_cells"],
            int(self.params["cil_mag"]), int(self.params["coa_mag"] * 16)
        )
        self.finished = False

    def confirm_int_steps_empty(self):
        if len(self.int_step_buffer) != 0:
            raise Exception("int steps buffer is not empty")

    def save_int_step(self, cell_ix, data):
        if self.current_cell_ix is not None:
            assert(self.current_cell_ix == cell_ix)
        else:
            self.current_cell_ix = cell_ix
        validate_data(self.data_labels, data)
        self.int_step_buffer.append(dict(copy.deepcopy(data)))
        if len(self.int_step_buffer) == self.num_int_steps:
            self.__save_cell_tstep(cell_ix)

    def __save_cell_tstep(self, cell_ix):
        validate_buffer(self.num_int_steps, self.int_step_buffer,
                        "INT_STEP_BUFFER")
        self.tstep_buffer[cell_ix] = copy.deepcopy(self.int_step_buffer)
        self.int_step_buffer = []
        self.current_cell_ix = None
        if np.all([len(self.tstep_buffer[ix]) == self.num_int_steps for ix in range(
                self.num_cells)]):
            self.__save_tstep()

    def __save_tstep(self):
        if not self.finished:
            validate_buffer(self.num_cells, self.tstep_buffer, "TSTEP_BUFFER")
            self.main_buffer["tsteps"].append(copy.deepcopy(self.tstep_buffer))
            self.tstep_buffer = [list() for _ in range(self.num_cells)]
            if len(self.main_buffer["tsteps"]) == self.num_tsteps:
                if os.path.exists(self.write_file_path):
                    os.remove(self.write_file_path)
                with open(self.write_file_path, "wb") as f:
                    f.write(orjson.dumps(self.main_buffer))
                self.finished = True
        else:
            raise Exception(
                "Supposed to be finished at {} tsteps.".format(self.num_tsteps))

    def save_header(self, data):
        data = make_params_list(data)
        validate_data(self.header_labels, data)
        self.main_buffer["header"] = dict(data)

    def close(self):
        if not self.finished:
            raise Exception("Expected {} tsteps, but only have {}.".format(
                self.num_tsteps, len(self.main_buffer["tsteps"])))
