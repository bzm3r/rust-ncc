from sim_data import SimulationData
from shared_sim_data import SharedSimData
from py_comp_data import PythonRustComparisonData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--exp", help="experiment json", type=str)
parser.add_argument("--output", help="output folder", type=str)
args = parser.parse_args()

exp_path = Path(args.exp)
exp_json = exp_path.parts[-1].replace(".json","") # Jank
with open(exp_path) as f:
    json_str = f.read()

exp_dict = orjson.loads(json_str)

out_dir = Path(args.output)

seeds, file_names = determine_file_names(exp_json, exp_dict)
for file_name in file_names:
    rust_dat = SimulationData()
    rust_dat.load_rust_dat(out_dir, file_name)
    rust_dat.tag = "rust"
    vec_ani_opts = get_vec_ani_opts(exp_dict)

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
