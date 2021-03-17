from sim_data import SimulationData
from sim_data import SharedSimData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson

run_experiments = True
exec_mode = "release"

root_dir = os.getcwd()
out_dir = os.path.join(root_dir, "output")
exp_jsons = ["two_cell_euler_ufine"]
poly_ls = ["-", ":"]

vec_ani_opts = []
sim_dats = []

for exp_json in exp_jsons:
    exec_path = os.path.join(root_dir, "target", exec_mode, "executor.exe")
    if run_experiments:
        build_out = subprocess.run(["cargo", "build"] + make_exec_mode_arg(
            exec_mode) + ["-p", "executor"])
        run_out = subprocess.run([exec_path] +
                                 ["-c", "./cfg.json"] +
                                 ["-e"] + exp_jsons)
        print(run_out)

    exp_path = os.path.join(root_dir, "experiments", "{}.json".format(exp_json))
    with open(exp_path) as f:
        json_str = f.read()

    exp_dict = orjson.loads(json_str)
    seeds, file_names = determine_file_names(exp_json, exp_dict)
    for (ix, file_name) in enumerate(file_names):
        sim_dat = SimulationData()
        sim_dat.load_rust_dat(out_dir, file_name)
