from paint_opts import AniOpts

def make_exec_mode_arg(exec_mode):
    if exec_mode == "debug":
        return []
    elif exec_mode == "release":
        return ["--release"]
    else:
        raise Exception("Unknown exec mode: {}".format(exec_mode))


def determine_file_names(exp_json, exp_dict):
    seeds = exp_dict["seeds"]

    file_names = []
    for seed in seeds:
        file_names.append("{}_seed={}".format(exp_json, seed))

    return seeds, file_names


def get_else(arg_dict, key, default):
    if key in arg_dict.keys():
        return arg_dict[key]
    else:
        return default


def get_vec_ani_opts(exp_dict):
    raw_vec_ani_opts = exp_dict["ani_opts"]
    vec_ani_opts = []
    for raw_ani_opts in raw_vec_ani_opts:
        rgtp_scale = get_else(raw_ani_opts, "rgtp_scale", 50)
        label_verts = get_else(raw_ani_opts, "label_verts", False)
        print("label_verts: ", label_verts)
        label_cells = get_else(raw_ani_opts, "label_cells", False)
        show_trails = get_else(raw_ani_opts, "show_trails", False)
        follow_group = get_else(raw_ani_opts, "follow_group", False)
        vec_ani_opts.append(AniOpts(rgtp_scale, label_verts, label_cells,
                                    show_trails, follow_group))
    return vec_ani_opts
