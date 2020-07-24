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

# circ_ixs = np.take(np.arange(16), np.arange(17), mode='wrap')
# fig, ax = plt.subplots()
#
# for i in range(len(vcs_per_vert_per_tstep)):
#     ax.cla()
#     ax.set_aspect('equal')
#     ax.set_xlim([-300, 300])
#     ax.set_ylim([-300, 300])
#     ax.plot(vcs_per_vert_per_tstep[i,circ_ixs,0], vcs_per_vert_per_tstep[i,circ_ixs,1], color='k')
#     ax.set_title("frame {}".format(i))
#     # Note that using time.sleep does *not* work here!
#     plt.pause(0.1)
