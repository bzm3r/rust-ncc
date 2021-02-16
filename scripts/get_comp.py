import hardio as hio
import os
import orjson
import copy

HEADER_LABELS = ["num_tsteps", "num_cells", "num_int_steps", "eta",
                 "f", "l", "t", "l3d", "kgtp", "close_zero_at",
                 "close_one_at", "cil_mag", "coa_los_penalty",
                 "coa_halfmax_dist",
                 "coa_mag", "coa_distrib_exp", "vertex_eta", "cell_r",
                 "rest_edge_len", "rest_area", "stiffness_edge",
                 "const_protrusive", "const_retractive", "stiffness_cyto",
                 "k_mem_on_vertex", "k_mem_off", "diffusion_rgtp", "init_rac",
                 "init_rho", "halfmax_vertex_rgtp_act",
                 "halfmax_vertex_rgtp_conc",
                 "kgtp_rac", "kgtp_rac_auto", "kdgtp_rac", "kdgtp_rho_on_rac",
                 "halfmax_tension_inhib", "tension_inhib", "kgtp_rho",
                 "kgtp_rho_auto", "kdgtp_rho", "kdgtp_rac_on_rho",
                 "randomization", "rand_avg_t", "rand_std_t", "rand_mag",
                 "num_rand_vs"]
BASICS = ["tstep", "int_step", "poly", "rac_acts", "rac_inacts", "rho_acts",
          "rho_inacts", "sum_forces"]
GEOMETRY = ["uivs"]
RAC_RATES = ["kgtps_rac", "kdgtps_rac"]
RHO_RATES = ["kgtps_rho", "kdgtps_rho"]
FORCES = ["rgtp_forces", "edge_forces", "cyto_forces"]
CALC_KGTPS_RAC = ["x_cils", "x_coas"]
OTHERS = ["edge_strains", "uevs", "poly_area", "coa_updates", "cil_updates"]

DATA_LABELS = BASICS + GEOMETRY + RAC_RATES + RHO_RATES + FORCES + \
              CALC_KGTPS_RAC + OTHERS


def read_save_file(out_dir, dat_file_path):
    fp = os.path.join(out_dir, dat_file_path)
    with open(fp, "rb") as f:
        out = orjson.loads(f.read())
    return out


# Input data has the form: data_per_int_step_per_cell_per_tstep
def get_data_per_c_per_s(data):
    num_int_steps = data["header"]["num_int_steps"]
    num_cells = data["header"]["num_cells"]
    tstep_data = data["tpoints"]
    num_tpoints = len(tstep_data)
    data_per_c_per_s = []
    for tix in range(num_tpoints):
        for ix in range(num_int_steps):
            data_per_c_per_s\
                .append(
                    [copy.deepcopy(
                        tstep_data[tix]["cells"][ci]["int_steps"][ix]
                    )
                     for ci in range(num_cells)]
                )
    return data_per_c_per_s
