from sim_data import SimulationData
from shared_sim_data import SharedSimData
from py_comp_data import PythonRustComparisonData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson

exp_name = "perftest_16_cell"
seeds = [7]
file_names = generate_file_names(seeds, exp_name)
root_dir = os.getcwd()
out_dir = os.path.join(root_dir, "output")
for file_name in file_names:
    rust_dat = SimulationData()
    rust_dat.load_rust_dat(out_dir, file_name)
    rust_dat.tag = "rust"
    vec_ani_opts = create_default_ani_opts()

    rust_dat.animate(vec_ani_opts, "rgtps")
    # rust_dat.animate(vec_ani_opts, "x_cils")
    # # rust_dat.animate(vec_ani_opts, "x_cals")
    # rust_dat.animate(vec_ani_opts, "kgtps_rho")
    # rust_dat.animate(vec_ani_opts, "kgtps_rac")
    # rust_dat.animate(vec_ani_opts, "rgtp_forces")
    # rust_dat.animate(vec_ani_opts, "x_coas")
    # rust_dat.animate(vec_ani_opts, "kdgtps_rac")
    # rust_dat.animate(vec_ani_opts, "kdgtps_rho")
    # rust_dat.animate(vec_ani_opts, "x_tens")
