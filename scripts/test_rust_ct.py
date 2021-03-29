from sim_data import SimulationData
from shared_sim_data import SharedSimData
from py_comp_data import PythonRustComparisonData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson
import copy

run_experiments = True
exec_mode = "release"
root_dir = os.getcwd()
char_ts = [2.0, 1.0, 0.5, 0.25, 0.1] #, 0.01, 0.005, 0.001]
exp_jsons = ["2_cell_ct_char_t_{}".format(x) for x in char_ts]
sim_dats = []
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
    for file_name in file_names:
        rust_dat = SimulationData()
        rust_dat.load_rust_dat(out_dir, file_name)
        rust_dat.tag = "rust"
        vec_ani_opts = get_vec_ani_opts(exp_dict)

        rust_dat.animate(vec_ani_opts, "rgtps")
        sim_dats.append(rust_dat)
        # rust_dat.animate(vec_ani_opts, "x_cils")
        # # rust_dat.animate(vec_ani_opts, "x_cals")
        # rust_dat.animate(vec_ani_opts, "kgtps_rho")
        # rust_dat.animate(vec_ani_opts, "kgtps_rac")
        # rust_dat.animate(vec_ani_opts, "rgtp_forces")
        # rust_dat.animate(vec_ani_opts, "x_coas")
        # rust_dat.animate(vec_ani_opts, "kdgtps_rac")
        # rust_dat.animate(vec_ani_opts, "kdgtps_rho")
        # rust_dat.animate(vec_ani_opts, "x_tens")

c0v0 = [sim_dat.poly_per_c_per_s[-1, 0, 0] for sim_dat in sim_dats]
c1v8 = [sim_dat.poly_per_c_per_s[-1, 1, 8] for sim_dat in sim_dats]
dists = np.array([np.linalg.norm(x - y) for x, y in zip(c0v0, c1v8)])
import matplotlib.pyplot as plt
plt.xlabel("dt")
plt.ylabel("final distance between vertices")
plt.semilogx(char_ts, dists, marker=".")
# out_dir = os.path.join(root_dir, "output")
# shared_sim_dat = SharedSimData(out_dir, sim_dats, poly_ls,
#                                "ct_combined")
# shared_sim_dat.animate(all_vec_ani_opts, "rgtps")