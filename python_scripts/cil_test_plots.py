import matplotlib.pyplot as plt
import numpy as np
import cbor2


def lookup_tstep_ix(tstep):
    return int(np.floor(tstep / frequency))


def v2ds_to_numpy(v2ds):
    vs = []
    for v2d in v2ds:
        vs.append([v2d['x'], v2d['y']])
    return np.array(vs)


def v2d_to_numpy(v2d):
    v = [v2d['x'], v2d['y']]
    return np.array(v)


def extract_v2ds_from_cell_states(state_key, dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['states']:
            dat_per_cell.append(v2ds_to_numpy(cell_rec[state_key][dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_v2ds_from_interactions(dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['interactions']:
            dat_per_cell.append(v2ds_to_numpy(cell_rec[dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_vector_to_close_points_from_interactions(state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['interactions']:
            x_close_points = cell_rec['x_close_points']
            dat_per_vertex = []
            for vertex in x_close_points:
                if vertex['empty']:
                    dat_per_vertex.append(np.array([np.nan, np.nan]))
                else:
                    vector_to = vertex['vector_to']
                    dat_per_vertex.append(v2d_to_numpy(vector_to))
            dat_per_cell.append(np.array(dat_per_vertex))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_scalars(state_key, dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['states']:
            dat_per_cell.append(np.array(cell_rec[state_key][dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


def extract_scalars_from_interactions(dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['interactions']:
            dat_per_cell.append(np.array(cell_rec[dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)


zs = ["true", "false"]
seeds = ["None"]
coas = ["24"]
cils = [60]
exp_names = ["pair"]

args = []
for z in zs:
    for coa in coas:
        for seed in seeds:
            for cil in cils:
                for exp_name in exp_names:
                    args.append([exp_name, cil, z, coa, seed])

for arg_set in args:
    exp_name, cil, z, coa, seed = arg_set
    file_path_template = "../output/{exp_name}_cil={cil}_Z={z}_cal=None_adh=None_coa={coa}_seed={seed}{extra}.{ext}"
    cbor_path = file_path_template.format(exp_name=exp_name, cil=cil, z=z, coa=coa, seed=seed, extra="", ext="cbor")
    graph_path = file_path_template.format(exp_name=exp_name, cil=cil, z=z, coa=coa, seed=seed, extra="{extra}", ext="png")
    snapshots = []
    with open(cbor_path, mode='rb') as sf:
        world_history = cbor2.load(sf)
        success = True
        while success:
            try:
                snapshots += cbor2.load(sf)
            except:
                success = False

    tsteps = [s["tstep"] for s in snapshots]
    state_recs = [s["cells"] for s in snapshots]
    frequency = world_history["snap_period"]

    poly_per_cell_per_tstep = extract_v2ds_from_cell_states('core', 'poly',
                                                            state_recs)
    centroids_per_cell_per_tstep = np.array(
        [[np.average(poly, axis=0) for poly in poly_per_cell] for poly_per_cell in
         poly_per_cell_per_tstep])
    rac_acts_per_cell_per_tstep = extract_scalars('core', 'rac_acts', state_recs)
    rho_acts_per_cell_per_tstep = extract_scalars('core', 'rho_acts', state_recs)

    adhs_per_cell_per_tstep = extract_v2ds_from_interactions('x_adhs', state_recs)
    rgtp_forces_per_cell_per_tstep = extract_v2ds_from_cell_states("mech",
                                                                   "rgtp_forces",
                                                                   state_recs)
    edge_forces_per_cell_per_tstep = extract_v2ds_from_cell_states("mech",
                                                                   "edge_forces",
                                                                   state_recs)
    cyto_forces_per_cell_per_tstep = extract_v2ds_from_cell_states("mech",
                                                                   "cyto_forces",
                                                                   state_recs)
    sum_non_adh_forces_per_cell_per_tstep = extract_v2ds_from_cell_states("mech",
                                                                          "sum_forces",
                                                                          state_recs)
    vec_to_cp_per_c_per_s = extract_vector_to_close_points_from_interactions(state_recs)
    mags_vec_to_cp_per_c_per_s = np.linalg.norm(vec_to_cp_per_c_per_s, axis=3)
    cil_per_c_per_s = extract_scalars_from_interactions('x_cils', state_recs)/60.0
    for vi in range(16):
        plt.plot(tsteps[:], cil_per_c_per_s[:, 0, vi], marker=".", ls="")
        plt.plot(tsteps[:], mags_vec_to_cp_per_c_per_s[:, 0, vi], marker=".", ls="")
        plt.title("seed={}, z={}, coa={}, ci={}, vi={}".format(seed, z, coa, 0, vi))
        plt.savefig(graph_path.format(extra="_ci={}_vi={}".format(0, vi)))
        plt.close()

