from sim_data import SimulationData
from sim_data import SharedSimData
from paint_opts import *
from utils import *
import os
import subprocess
import toml


root_dir = os.getcwd()
exp_tomls = ["py_comp_2_long", "py_comp_2_long_reduced_cil"]
for exp_toml in exp_tomls:
    exec_mode = "release"
    exec_path = os.path.join(root_dir, "target", exec_mode, "executor.exe")
    if True:
        build_out = subprocess.run(["cargo", "build"] + make_exec_mode_arg(
            exec_mode) + ["-p", "executor"])
        run_out = subprocess.run([exec_path] +
                                 ["-c", "./cfg.toml"] +
                                 ["-e"] + exp_tomls)
        print(run_out)

    exp_dict = toml.load(
        os.path.join(root_dir, "experiments", "{}.toml".format(exp_toml))
    )
    file_names = determine_file_names(exp_toml, exp_dict)
    out_dir = os.path.join(root_dir, "output")
    for file_name in file_names:
        rust_dat = SimulationData()
        rust_dat.load_rust_dat(out_dir, file_name)
        rust_dat.tag = "rust"
        py_dat = SimulationData()
        py_dat.load_py_dat(out_dir, file_name)
        py_dat.tag = "python"
        vec_ani_opts = get_vec_ani_opts(exp_dict)
        # rust_dat.animate(vec_ani_opts)
        # py_dat.animate(vec_ani_opts)
        comp_dat = SharedSimData(out_dir, [rust_dat, py_dat], ["-", ":"],
                                 file_name +
                                 "_rust_and_py")
        comp_dat.animate(vec_ani_opts)
