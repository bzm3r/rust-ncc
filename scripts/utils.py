from paint_opts import AniOpts


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


def get_else(toml_dict, key, default):
    if key in toml_dict.keys():
        return toml_dict[key]
    else:
        return default


def get_vec_ani_opts(exp_dict):
    vec_raw_ani_opts = get_else(exp_dict, "ani-opts", [])
    vec_ani_opts = []
    for raw_ani_opts in vec_raw_ani_opts:
        rgtp_scale = get_else(raw_ani_opts, "rgtp_scale", 50)
        label_verts = get_else(raw_ani_opts, "label_verts", False)
        print("label_verts: ", label_verts)
        label_cells = get_else(raw_ani_opts, "label_cells", False)
        show_trails = get_else(raw_ani_opts, "show_trails", False)
        follow_group = get_else(raw_ani_opts, "follow_group", False)
        vec_ani_opts.append(AniOpts(rgtp_scale, label_verts, label_cells,
                                    show_trails, follow_group))
    return vec_ani_opts
