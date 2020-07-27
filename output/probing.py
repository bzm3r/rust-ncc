import matplotlib.pyplot as plt
import fastavro
import json
import numpy as np
import copy

state_recs = []
with open('history_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('history_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            state_recs.append(r)

poly_per_tstep = []
for rec in state_recs:
    poly = []
    for xy in rec['cells'][0]['state']['vertex_coords']:
        poly.append([xy['x'], xy['y']])
    poly_per_tstep.append(copy.deepcopy(poly))
poly_per_tstep = np.array(poly_per_tstep)

geom_recs = []
with open('geom_hist_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('geom_hist_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            geom_recs.append(r)

mech_recs = []
with open('mech_hist_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('mech_hist_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            mech_recs.append(r)

uivs_per_tstep = []
for poly, rec in zip(poly_per_tstep, geom_recs):
    uivs = []
    for vc, uiv in zip(poly, rec['state'][0]['unit_inward_vecs']):
        uivs.append([uiv['x'], uiv['y']])
    uivs_per_tstep.append(copy.deepcopy(uivs))
uivs_per_tstep = np.array(uivs_per_tstep)

uevs_per_tstep = []
for poly, rec in zip(poly_per_tstep, geom_recs):
    uevs = []
    for vc, uev in zip(poly, rec['state'][0]['unit_edge_vecs']):
        uevs.append([uev['x'], uev['y']])
    uevs_per_tstep.append(copy.deepcopy(uevs))
uevs_per_tstep = np.array(uevs_per_tstep)

rac_acts_arrows_per_tstep = []
rac_acts_per_tstep = []
for tstep, rec in enumerate(state_recs):
    rac_acts = []
    rac_acts_arrows = []
    for (vix, x) in enumerate(rec['cells'][0]['state']['rac_acts']):
        rac_acts.append(x)
        arrow_deltas = -1500 * x * uivs_per_tstep[tstep][vix]
        rac_acts_arrows.append([poly_per_tstep[tstep][vix][0] + uevs_per_tstep[tstep][vix][0]*0.1, poly_per_tstep[tstep][vix][1] + uevs_per_tstep[tstep][vix][0]*0.1, arrow_deltas[0], arrow_deltas[1]])
    rac_acts_arrows_per_tstep.append(copy.deepcopy(rac_acts_arrows))
    rac_acts_per_tstep.append(rac_acts)
rac_acts_per_tstep = np.array(rac_acts_per_tstep)
rac_acts_arrows_per_tstep = np.array(rac_acts_arrows_per_tstep)

rho_acts_arrows_per_tstep = []
rho_acts_per_tstep = []
for tstep, rec in enumerate(state_recs):
    rho_acts = []
    rho_acts_arrows = []
    for (vix, x) in enumerate(rec['cells'][0]['state']['rho_acts']):
        rho_acts.append(x)
        arrow_deltas = -1500 * x * uivs_per_tstep[tstep][vix]
        rho_acts_arrows.append([poly_per_tstep[tstep][vix][0] + uevs_per_tstep[tstep][vix][0]*0.1, poly_per_tstep[tstep][vix][1] + uevs_per_tstep[tstep][vix][0]*0.1, arrow_deltas[0], arrow_deltas[1]])
    rho_acts_arrows_per_tstep.append(copy.deepcopy(rho_acts_arrows))
    rho_acts_per_tstep.append(rho_acts)
rho_acts_per_tstep = np.array(rho_acts_per_tstep)
rho_acts_arrows_per_tstep = np.array(rho_acts_arrows_per_tstep)

cyto_forces_per_tstep = []
for tstep, rec in enumerate(mech_recs):
    cyto_forces = []
    for (vix, xy) in enumerate(rec['state'][0]['cyto_forces']):
        arrow_deltas = 0.05*np.array([xy['x'], xy['y']])
        cyto_forces.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
    cyto_forces_per_tstep.append(copy.deepcopy(cyto_forces))
cyto_forces_per_tstep = np.array(cyto_forces_per_tstep)

edge_forces_plus_per_tstep = []
for tstep, rec in enumerate(mech_recs):
    efs_plus = []
    for vix in range(16):
        xy = rec['state'][0]['edge_forces'][vix]
        arrow_deltas = 0.1*np.array([xy['x'], xy['y']])
        efs_plus.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
    edge_forces_plus_per_tstep.append(copy.deepcopy(efs_plus))
edge_forces_plus_per_tstep = np.array(edge_forces_plus_per_tstep)

edge_forces_minus_per_tstep = []
for tstep, rec in enumerate(mech_recs):
    efs_minus = []
    for vix in range(16):
        xy = rec['state'][0]['edge_forces'][(vix - 1)%16]
        arrow_deltas = -0.1*np.array([xy['x'], xy['y']])
        efs_minus.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
    edge_forces_minus_per_tstep.append(copy.deepcopy(efs_minus))
edge_forces_minus_per_tstep = np.array(edge_forces_minus_per_tstep)

edge_forces_per_tstep = np.append(poly_per_tstep, edge_forces_plus_per_tstep[:,:,2:4] - edge_forces_minus_per_tstep[:,:,2:4], axis=2)

edge_strains_per_tstep = []
for tstep, rec in enumerate(mech_recs):
    edge_strains = []
    for vix in range(16):
        edge_strains.append(rec['state'][0]['edge_strains'][vix])
    edge_strains_per_tstep.append(copy.deepcopy(edge_strains))
edge_strains_per_tstep = np.array(edge_strains_per_tstep)

rgtp_forces_per_tstep = []
for tstep, rec in enumerate(mech_recs):
    rgtp_forces = []
    for (vix, xy) in enumerate(rec['state'][0]['rgtp_forces']):
        arrow_deltas = 0.05*np.array([xy['x'], xy['y']])
        rgtp_forces.append([poly_per_tstep[tstep][vix][0], poly_per_tstep[tstep][vix][1], arrow_deltas[0], arrow_deltas[1]])
    rgtp_forces_per_tstep.append(copy.deepcopy(rgtp_forces))
rgtp_forces_per_tstep = np.array(rgtp_forces_per_tstep)

circ_vixs = np.take(np.arange(16), np.arange(17), mode='wrap')
def paint(delta):
    global fig
    global ax
    global tstep
    global num_tsteps
    ax.cla()
    ax.set_aspect('equal')
    ax.set_xlim([20.0 - 150.0, 20.0 + 150.0])
    ax.set_ylim([20.0 - 150.0, 20.0 + 150.0])
    for vix in range(16):
        print(edge_strains_per_tstep[tstep])
        if False:#edge_strains_per_tstep[tstep][vix] > 0.0:
            strain_ls = ":"
        else:
            strain_ls = "-"
        ax.plot([poly_per_tstep[tstep, vix, 0], poly_per_tstep[tstep, (vix + 1)%16, 0]], [poly_per_tstep[tstep, vix, 1], poly_per_tstep[tstep, (vix + 1)%16, 1]], color='k', ls=strain_ls)
    for rac_act in rac_acts_arrows_per_tstep[tstep]:
        ax.arrow(rac_act[0], rac_act[1], rac_act[2], rac_act[3], color="b", length_includes_head=True, head_width=0.0)
    for rho_act in rho_acts_arrows_per_tstep[tstep]:
        ax.arrow(rho_act[0], rho_act[1], rho_act[2], rho_act[3], color="r", length_includes_head=True, head_width=0.0)
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
    ax.set_title("frame {}".format(tstep))
    tstep = (tstep + delta) % num_tsteps
    plt.show()


def on_press(event):
    global fig
    if event.key == 'x':
        paint(1)
    elif event.key == 'z':
        paint(-1)
    elif event.key == 'q':
        paint(-100)
    elif event.key == 'w':
        paint(100)
    fig.canvas.draw()


num_tsteps = poly_per_tstep.shape[0]
tstep = 0
fig, ax = plt.subplots()
fig.canvas.mpl_connect('key_press_event', on_press)
