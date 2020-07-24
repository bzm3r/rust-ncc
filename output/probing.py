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

mech_recs = []
with open('mech_hist_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('mech_hist_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            mech_recs.append(r)

geom_recs = []
with open('geom_hist_schema.avsc') as sf:
    fastavro.parse_schema(json.load(sf))
    with open('geom_hist_dat.avro', 'rb') as df:
        for r in fastavro.reader(df):
            geom_recs.append(r)

ef_per_vert_per_tstep = []
for poly, rec in zip(poly_per_tstep, mech_recs):
    efs = []
    arrow_width = 1
    for vc, ef in zip(poly, rec['state'][0]['edge_forces']):
        efs.append([vc[0], vc[1], 10*ef['x'], 10*ef['y']])
    ef_per_vert_per_tstep.append(copy.deepcopy(efs))
ef_per_vert_per_tstep = np.array(ef_per_vert_per_tstep)

circ_ixs = np.take(np.arange(16), np.arange(17), mode='wrap')
fig, ax = plt.subplots()

for i in range(len(poly_per_tstep)):
    ax.cla()
    ax.set_aspect('equal')
    ax.set_xlim([-75, 75])
    ax.set_ylim([-75, 75])
    ax.plot(poly_per_tstep[i,circ_ixs,0], poly_per_tstep[i,circ_ixs,1], color='k')
    for ef in ef_per_vert_per_tstep[i]:
        print(ef)
        ax.arrow(ef[0], ef[1], ef[2], ef[3], color="g", length_includes_head=True, head_width=1)
    ax.set_title("frame {}".format(i))
    # Note that using time.sleep does *not* work here!
    plt.pause(2.5)
