from sim_data import SimulationData
from paint_opts import *
import os
import subprocess
import toml


def get_optional(toml_dict, key, default):
    if key in toml_dict.keys():
        return toml_dict[key]
    else:
        return default


def get_seeds(toml_dict):
    raw_seeds = get_optional(toml_dict, "seeds", [])
    seeds = []
    for rs in raw_seeds:
        seeds.append(get_optional(rs, "seed", None))
    return seeds


def make_exec_mode_arg(exec_mode):
    if exec_mode == "debug":
        return []
    elif exec_mode == "release":
        return ["--release"]
    else:
        raise Exception("Unknown exec mode: {}".format(exec_mode))


def determine_file_names(exp_toml, exp_dict):
    seeds = get_seeds(exp_dict)

    file_names = []
    for seed in seeds:
        file_names.append("{}_seed={}".format(exp_toml, seed))

    return file_names


def get_else_default(toml_dict, key, default):
    if key in toml_dict.keys():
        return toml_dict[key]
    else:
        return default


def get_vec_ani_opts(exp_dict):
    vec_raw_ani_opts = get_else_default(exp_dict, "ani-opts", [])
    vec_ani_opts = []
    for raw_ani_opts in vec_raw_ani_opts:
        rgtp_scale = get_else_default(raw_ani_opts, "rgtp_scale", 50)
        label_verts = get_else_default(raw_ani_opts, "label_verts", False)
        label_cells = get_else_default(raw_ani_opts, "label_cells", False)
        show_trails = get_else_default(raw_ani_opts, "show_trails", False)
        follow_group = get_else_default(raw_ani_opts, "follow_group", False)
        vec_ani_opts.append(AniOpts(rgtp_scale, label_verts, label_cells,
                                    show_trails, follow_group))
    return vec_ani_opts


root_dir = "B:\\rust-ncc\\"
exp_tomls = ["pair_2_short"]  # ["test", "py_comp_1", "py_comp_2"]
for exp_toml in exp_tomls:
    exec_mode = "debug"
    exec_path = os.path.join(root_dir, "target", exec_mode, "executor.exe")
    if False:
        build_out = subprocess.run(["cargo", "build"] + make_exec_mode_arg(
            exec_mode) + ["-p", "executor"])
        run_out = subprocess.run([exec_path] + exp_tomls)
        print(run_out)

    exp_dict = toml.load(
        os.path.join(root_dir, "experiments", "{}.toml".format(exp_toml))
    )
    file_names = determine_file_names(exp_toml, exp_dict)
    out_dir = os.path.join(root_dir, "output")
    for file_name in file_names:
        sim_dat = SimulationData(out_dir, file_name)
        vec_ani_opts = get_vec_ani_opts(exp_dict)
        sim_dat.animate(vec_ani_opts)
