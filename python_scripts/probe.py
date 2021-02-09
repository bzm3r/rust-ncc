import matplotlib.pyplot as plt
import json
import numpy as np
import json
import cbor2


output = None

file_name = "../output/separated_pair_cil=60_cal=None_adh=10_coa=24_seed=9757.cbor"

snapshots = []
with open(file_name, mode='rb') as sf:
    world_history = cbor2.load(sf)
    success = True
    while success:
        try:
            snapshots += cbor2.load(sf)
        except:
            success = False

tsteps = [s["tstep"] for s in snapshots]
state_recs = [s["cells"] for s in snapshots]
frequency = world_history["snap_freq"]


def lookup_tstep_ix(tstep):
    return int(np.floor(tstep/frequency))

def p2ds_to_numpy(p2ds):
    vs = []
    for p2d in p2ds:
        vs.append([p2d['x'], p2d['y']])
    return np.array(vs)


def extract_p2ds_from_cell_states(state_key, dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['states']:
            dat_per_cell.append(p2ds_to_numpy(cell_rec[state_key][dat_key]))
        dat_per_cell_per_tstep.append(np.array(dat_per_cell))
    return np.array(dat_per_cell_per_tstep)

def extract_p2ds_from_interactions(dat_key, state_recs):
    dat_per_cell_per_tstep = []
    for rec in state_recs:
        dat_per_cell = []
        for cell_rec in rec['interactions']:
            dat_per_cell.append(p2ds_to_numpy(cell_rec[dat_key]))
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


poly_per_cell_per_tstep = extract_p2ds_from_cell_states('core', 'poly', state_recs)
centroids_per_cell_per_tstep = np.array(
    [[np.average(poly, axis=0) for poly in poly_per_cell] for poly_per_cell in
     poly_per_cell_per_tstep])
uivs_per_cell_per_tstep = extract_p2ds_from_cell_states('geom', 'unit_inward_vecs',
                                                        state_recs)
uovs_per_cell_per_tstep = -1 * uivs_per_cell_per_tstep
rac_acts_per_cell_per_tstep = extract_scalars('core', 'rac_acts', state_recs)
rac_act_arrows_per_cell_per_tstep = 50 * rac_acts_per_cell_per_tstep[:, :, :,
                                         np.newaxis] * uovs_per_cell_per_tstep
rho_acts_per_cell_per_tstep = extract_scalars('core', 'rho_acts', state_recs)
rho_act_arrows_per_cell_per_tstep = 50 * rho_acts_per_cell_per_tstep[:, :, :,
                                         np.newaxis] * uivs_per_cell_per_tstep

adhs_per_cell_per_tstep = 5*extract_p2ds_from_interactions('x_adhs', state_recs)

#
# rho_acts_arrows_per_tstep = []
# rho_acts_per_tstep = []
# for tstep, rec in enumerate(state_recs):
#     rho_acts = []
#     rho_acts_arrows = []
#     for (vix, x) in enumerate(rec['cells'][0]['state']['rho_acts']):
#         rho_acts.append(x)
#         arrow_deltas = -150 * x * uivs_per_tstep[tstep][vix]
#         rho_acts_arrows.append([poly_per_tstep[tstep][vix][0] + uevs_per_tstep[tstep][vix][0]*0.1, poly_per_tstep[tstep][vix][1] + uevs_per_tstep[tstep][vix][0]*0.1, arrow_deltas[0], arrow_deltas[1]])
#     rho_acts_arrows_per_tstep.append(copy.deepcopy(rho_acts_arrows))
#     rho_acts_per_tstep.append(rho_acts)
# rho_acts_per_tstep = np.array(rho_acts_per_tstep)
# rho_acts_arrows_per_tstep = np.array(rho_acts_arrows_per_tstep)
#
# cyto_forces_per_tstep = []
# for tstep, rec in enumerate(mech_recs):
#     cyto_forces = []
#     for (vix, xy) in enumerate(rec['state'][0]['cyto_forces']):
#         arrow_deltas = 0.05*np.array([xy['x'], xy['y']])
#         cyto_forces.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
#     cyto_forces_per_tstep.append(copy.deepcopy(cyto_forces))
# cyto_forces_per_tstep = np.array(cyto_forces_per_tstep)
#
# edge_forces_plus_per_tstep = []
# for tstep, rec in enumerate(mech_recs):
#     efs_plus = []
#     for vix in range(16):
#         xy = rec['state'][0]['edge_forces'][vix]
#         arrow_deltas = 0.1*np.array([xy['x'], xy['y']])
#         efs_plus.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
#     edge_forces_plus_per_tstep.append(copy.deepcopy(efs_plus))
# edge_forces_plus_per_tstep = np.array(edge_forces_plus_per_tstep)
#
# edge_forces_minus_per_tstep = []
# for tstep, rec in enumerate(mech_recs):
#     efs_minus = []
#     for vix in range(16):
#         xy = rec['state'][0]['edge_forces'][(vix - 1)%16]
#         arrow_deltas = -0.1*np.array([xy['x'], xy['y']])
#         efs_minus.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
#     edge_forces_minus_per_tstep.append(copy.deepcopy(efs_minus))
# edge_forces_minus_per_tstep = np.array(edge_forces_minus_per_tstep)
#
# edge_forces_per_tstep = np.append(poly_per_tstep, edge_forces_plus_per_tstep[:,:,2:4] - edge_forces_minus_per_tstep[:,:,2:4], axis=2)
#
# edge_strains_per_tstep = []
# for tstep, rec in enumerate(mech_recs):
#     edge_strains = []
#     for vix in range(16):
#         edge_strains.append(rec['state'][0]['edge_strains'][vix])
#     edge_strains_per_tstep.append(copy.deepcopy(edge_strains))
# edge_strains_per_tstep = np.array(edge_strains_per_tstep)
#
# rgtp_forces_per_tstep = []
# for tstep, rec in enumerate(mech_recs):
#     rgtp_forces = []
#     for (vix, xy) in enumerate(rec['state'][0]['rgtp_forces']):
#         arrow_deltas = 0.05*np.array([xy['x'], xy['y']])
#         rgtp_forces.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
#     rgtp_forces_per_tstep.append(copy.deepcopy(rgtp_forces))
# rgtp_forces_per_tstep = np.array(rgtp_forces_per_tstep)

circ_vixs = np.take(np.arange(16), np.arange(17), mode='wrap')


def paint(delta):
    global fig
    global ax
    global tstep_ix
    global num_tsteps

    old_xlim = ax.get_xlim()
    old_ylim = ax.get_ylim()
    ax.cla()
    ax.set_xlim(old_xlim)
    ax.set_ylim(old_ylim)

    for (ci, poly) in enumerate(poly_per_cell_per_tstep[tstep_ix]):
        if ci == 0:
            poly_color = "k"
        else:
            poly_color = "g"

        for vix in range(16):
            ax.plot([poly[vix, 0], poly[(vix + 1) % 16, 0]],
                    [poly[vix, 1], poly[(vix + 1) % 16, 1]],
                    color=poly_color, marker=".")
            ax.annotate(str(vix), (poly[vix, 0], poly[vix, 1]))

    for poly, rac_act_arrows in zip(
            poly_per_cell_per_tstep[tstep_ix],
            rac_act_arrows_per_cell_per_tstep[tstep_ix]
    ):
        for p, rac_arrow in zip(poly, rac_act_arrows):
            ax.arrow(p[0], p[1], 1*rac_arrow[0], 1*rac_arrow[1], color="b",
                     length_includes_head=True, head_width=0.0)

    for poly, rho_act_arrows in zip(poly_per_cell_per_tstep[tstep_ix],
                                    rho_act_arrows_per_cell_per_tstep[tstep_ix]):
        for p, rho_arrow in zip(poly, rho_act_arrows):
            ax.arrow(p[0], p[1], 1*rho_arrow[0], 1*rho_arrow[1], color="r",
                     length_includes_head=True, head_width=0.0)

    for poly_ix, poly, adhs in zip(np.arange(0, len(poly_per_cell_per_tstep[0])), poly_per_cell_per_tstep[tstep_ix], adhs_per_cell_per_tstep[tstep_ix]):
        if poly_ix == 0:
            adh_arrow_color = "magenta"
        else:
            adh_arrow_color = "cyan"
        for p, adh in zip(poly, adhs):
            ax.arrow(p[0], p[1], adh[0], adh[1], color=adh_arrow_color,
             length_includes_head=True, head_width=1.0)

    # for rac_act in rac_acts_arrows_per_tstep[tstep]:
    #     ax.arrow(rac_act[0], rac_act[1], rac_act[2], rac_act[3], color="b", length_includes_head=True, head_width=0.0)
    # for rho_act in rho_acts_arrows_per_tstep[tstep]:
    #     ax.arrow(rho_act[0], rho_act[1], rho_act[2], rho_act[3], color="r", length_includes_head=True, head_width=0.0)
    # for cyto_force in cyto_forces_per_tstep[tstep]:
    #     ax.arrow(cyto_force[0], cyto_force[1], cyto_force[2], cyto_force[3], color="cyan", length_includes_head=True, head_width=0.5)
    # for i, ef_plus in enumerate(edge_forces_plus_per_tstep[tstep]):
    #     if True:
    #         ax.plot([ef_plus[0]], [ef_plus[1]], marker='o', color='b')
    #         ax.arrow(ef_plus[0], ef_plus[1], ef_plus[2], ef_plus[3], color="b", length_includes_head=True, head_width=0.5)
    # for i, ef_minus in enumerate(edge_forces_minus_per_tstep[tstep]):
    #     if True:
    #         ax.plot([ef_minus[0]], [ef_minus[1]], marker='o', color='r')
    #         ax.arrow(ef_minus[0], ef_minus[1], ef_minus[2], ef_minus[3], color="r", length_includes_head=True, head_width=0.5)
    # for edge_force in edge_forces_per_tstep[tstep]:
    #     ax.arrow(edge_force[0], edge_force[1], edge_force[2], edge_force[3], color="g", length_includes_head=True, head_width=0.5)
    # for rgtp_force in rgtp_forces_per_tstep[tstep]:
    #     ax.arrow(rgtp_force[0], rgtp_force[1], rgtp_force[2], rgtp_force[3], color="b", length_includes_head=True, head_width=0.5)
    ax.set_title("frame {}".format(tsteps[tstep_ix]))
    tstep_ix = (tstep_ix + delta) % len(tsteps)
    plt.show()


def on_press(event):
    global fig
    global DEFAULT_XLIM
    global DEFAULT_YLIM
    if event.key == 'x':
        paint(1)
    elif event.key == 'z':
        paint(-1)
    if event.key == 'c':
        paint(-5)
    elif event.key == 'v':
        paint(5)
    elif event.key == 'n':
        paint(-10)
    elif event.key == 'm':
        paint(10)
    elif event.key == 'r':
        ax.set_aspect('equal')
        ax.set_xlim(DEFAULT_XLIM)
        ax.set_ylim(DEFAULT_YLIM)
    fig.canvas.draw()


num_tsteps = poly_per_cell_per_tstep.shape[0]
tstep_ix = 0
fig, ax = plt.subplots()
DEFAULT_XLIM = [-100, 200]
DEFAULT_YLIM = [-100, 200]
ax.set_aspect('equal')
ax.set_xlim(DEFAULT_XLIM)
ax.set_ylim(DEFAULT_YLIM)
fig.canvas.mpl_connect('key_press_event', on_press)
