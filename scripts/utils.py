from paint_opts import AniOpts
import numpy as np


def make_exec_mode_arg(exec_mode):
    if exec_mode == "debug":
        return []
    elif exec_mode == "release":
        return ["--release"]
    else:
        raise Exception("Unknown exec mode: {}".format(exec_mode))


def determine_file_names(exp_json, exp_dict):
    seeds = exp_dict["seeds"]

    file_names = generate_file_names(seeds, exp_json)

    return seeds, file_names

def generate_file_names(seeds, exp_name):
    file_names = []
    for seed in seeds:
        file_names.append("{}_seed={}".format(exp_name, seed))

    return file_names


def get_else(arg_dict, key, default):
    if key in arg_dict.keys():
        return arg_dict[key]
    else:
        return default


def get_vec_ani_opts(exp_dict):
    raw_vec_ani_opts = exp_dict["ani_opts"]
    vec_ani_opts = []
    for raw_ani_opts in raw_vec_ani_opts:
        arrow_scale = get_else(raw_ani_opts, "arrow_scale", 50)
        label_verts = get_else(raw_ani_opts, "label_verts", False)
        print("label_verts: ", label_verts)
        label_cells = get_else(raw_ani_opts, "label_cells", False)
        show_trails = get_else(raw_ani_opts, "show_trails", False)
        follow_group = get_else(raw_ani_opts, "follow_group", False)
        vec_ani_opts.append(AniOpts(arrow_scale, label_verts, label_cells,
                                    show_trails, follow_group))
    return vec_ani_opts


def create_default_ani_opts():
    vec_ani_opts = []

    arrow_scale = 1
    label_verts = False
    label_cells = False
    show_trails = False
    follow_group = False
    vec_ani_opts.append(AniOpts(arrow_scale, label_verts, label_cells,
                                show_trails, follow_group))
    return vec_ani_opts


def find_common_ts(ixs_ts_per_sim):
    if len(ixs_ts_per_sim) == 0:
        return []
    elif len(ixs_ts_per_sim) == 1:
        return ([t for (ix, t) in ixs_ts_per_sim[0]],
                [[ix for (ix, t) in ixs_ts] for ixs_ts in ixs_ts_per_sim])
    else:
        common = []
        ixs_per_list = [list() for _ in range(len(ixs_ts_per_sim))]
        num_commons = 0
        ixs_ts_per_sim = sorted(ixs_ts_per_sim, key=lambda x: len(x))
        xs = ixs_ts_per_sim[0]
        num_xs = len(xs)
        min_ixs = [0 for _ in ixs_ts_per_sim[1:]]
        for (x_snap_ix, x) in xs:
            ix_per_list = [x_snap_ix]
            if num_commons < num_xs:
                inside_all = True
                for (list_ix, ys) in enumerate(ixs_ts_per_sim[1:]):
                    if inside_all:
                        check = False
                        for (ix_y, (y_snap_ix, y)) in enumerate(ys):
                            if abs((x - y)) < 1e-4:
                                check = True
                                ix_per_list.append(y_snap_ix)
                                break
                        inside_all = check and inside_all
                    else:
                        break
                if inside_all:
                    common.append(x)
                    if len(ix_per_list) != len(ixs_per_list):
                        raise Exception(
                            "len(ix_per_list) = {} != {} = len(ixs_per_list)"
                                .format(len(ix_per_list), len(ixs_per_list)))
                    for (ix, ixs) in zip(ix_per_list, ixs_per_list):
                        ixs.append(ix)

                    num_commons += 1

        return common, ixs_per_list


def sanitize_tpoints(tpoints):
    tpoints = np.array(tpoints)
    deltas = tpoints[1:] - tpoints[:-1]
    max_delta = np.max(deltas)
    is_max = np.flatnonzero(np.abs(deltas - max_delta) < 1e-4) + 1
    sanitized_tpoints = [tpoints[0]] + tpoints[is_max].tolist()
    snap_ixs = [0] + is_max.tolist()
    return list(zip(snap_ixs, sanitized_tpoints))


