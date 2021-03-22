from sim_data import SimulationData
from shared_sim_data import SharedSimData
from py_comp_data import PythonRustComparisonData
from paint_opts import *
from utils import *
import os
import subprocess
import orjson

run_experiments = True
exec_mode = "release"
root_dir = os.getcwd()
exp_jsons = ["py_comp_2"]
for exp_json in exp_jsons:
    exec_path = os.path.join(root_dir, "target", exec_mode, "executor")
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
    out_dir = os.path.join(root_dir, "output")
    for file_name in file_names:
        rust_dat = SimulationData()
        rust_dat.load_rust_dat(out_dir, file_name)
        py_dat = SimulationData()
        py_dat.load_py_dat(out_dir, file_name)
        vec_ani_opts = get_vec_ani_opts(exp_dict)

        rust_dat.animate(vec_ani_opts, "rgtps")
        # rust_dat.animate(vec_ani_opts, "rho_acts")
        # rust_dat.animate(vec_ani_opts, "x_cils")
        # rust_dat.animate(vec_ani_opts, "kgtps_rho")
        # rust_dat.animate(vec_ani_opts, "kgtps_rac")
        # rust_dat.animate(vec_ani_opts, "rgtp_forces")
        # rust_dat.animate(vec_ani_opts, "x_coas")
        # rust_dat.animate(vec_ani_opts, "kdgtps_rac")
        # rust_dat.animate(vec_ani_opts, "kdgtps_rho")
        # rust_dat.animate(vec_ani_opts, "rho_act_net_fluxes")
        # rust_dat.animate(vec_ani_opts, "rho_inacts")
        #
        py_dat.animate(vec_ani_opts, "rgtps")
        # py_dat.animate(vec_ani_opts, "rho_acts")
        # py_dat.animate(vec_ani_opts, "x_cils")
        # py_dat.animate(vec_ani_opts, "kgtps_rho")
        # py_dat.animate(vec_ani_opts, "kgtps_rac")
        # py_dat.animate(vec_ani_opts, "rgtp_forces")
        # py_dat.animate(vec_ani_opts, "x_coas")
        # py_dat.animate(vec_ani_opts, "kdgtps_rac")
        # py_dat.animate(vec_ani_opts, "kdgtps_rho")
        # py_dat.animate(vec_ani_opts, "rho_act_net_fluxes")
        # py_dat.animate(vec_ani_opts, "rho_inacts")

        comp_dat = PythonRustComparisonData(out_dir, py_dat, rust_dat,
                                            [":", "-"], file_name +
                                            "_rust_and_py")
        comp_dat.plot(["kdgtps_rac", "kdgtps_rho", "rho_acts", "kgtps_rho",
                       "kdgtps_rho"])
        comp_dat.animate(vec_ani_opts, "rgtps")
        # comp_dat.animate(vec_ani_opts, "rho_acts")
        # comp_dat.animate(vec_ani_opts, "kdgtps_rho")
        # comp_dat.animate(vec_ani_opts, "kgtps_rho")
        # comp_dat.animate(vec_ani_opts, "kdgtps_rac")
        # comp_dat.animate(vec_ani_opts, "kgtps_rac")
        # comp_dat.animate(vec_ani_opts, "rho_act_net_fluxes")
        # rust_dat.animate(vec_ani_opts, "rho_inacts")
