from sim_data import SimulationData
from shared_sim_data import SharedSimData
from py_comp_data import PythonRustComparisonData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson

run_experiments = False
exec_mode = "release"
root_dir = os.getcwd()
exp_jsons = ["two_cell_coa_12", "two_cell_coa_24"]
poly_ls = [":", "-"]
sim_dats = []
all_vec_ani_opts = []
for exp_json in exp_jsons:
    exec_path = os.path.join(root_dir, "target", exec_mode, "executor")
    if run_experiments:
        build_out = subprocess.run(["cargo", "build"] + make_exec_mode_arg(
            exec_mode) + ["-p", "executor"])
        run_out = subprocess.run([exec_path] +
                                 ["-c", "./cfg.json"] +
                                 ["-e"] + [exp_json])
        print(run_out)

    exp_path = os.path.join(root_dir, "experiments", "{}.json".format(exp_json))
    with open(exp_path) as f:
        json_str = f.read()

    exp_dict = orjson.loads(json_str)
    seeds, file_names = determine_file_names(exp_json, exp_dict)
    out_dir = os.path.join(root_dir, "output")
    for (ix, file_name) in enumerate(file_names):
        rust_dat = SimulationData()
        rust_dat.load_rust_dat(out_dir, file_name)
        rust_dat.tag = "rust"
        vec_ani_opts = get_vec_ani_opts(exp_dict)
        if ix == 0:
            all_vec_ani_opts = vec_ani_opts
            sim_dats.append(rust_dat)
        # rust_dat.animate(vec_ani_opts, "rgtps")

shared_sim_dat = SharedSimData(out_dir, sim_dats, poly_ls,
                               "combined_animation")
shared_sim_dat.animate(all_vec_ani_opts, "rgtps")
