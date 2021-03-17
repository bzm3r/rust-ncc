from sim_data import SimulationData
from sim_data import SharedSimData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson
import copy

run_experiments = False
exec_mode = "release"

root_dir = os.getcwd()
out_dir = os.path.join(root_dir, "output")
exp_jsons = ["two_cell_euler_coarse", "two_cell_euler_medium",
             "two_cell_euler_fine", "two_cell_euler_ufine"][2:]
poly_ls = [":", "-.", "--", "-"][2:]

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
        sim_dat.load_dat(out_dir, file_name)
        if ix == 0:
            vec_ani_opts = get_vec_ani_opts(exp_dict)
            sim_dats.append(copy.deepcopy(sim_dat))

comp_dat = SharedSimData(out_dir, sim_dats, poly_ls, "euler_int_step_test_fine")
comp_dat.animate(vec_ani_opts)
