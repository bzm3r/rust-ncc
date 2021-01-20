import os
import numpy as np
import cbor2
import parse
import copy

test_file = "../output/out_py.dat"

delimiter = "++++++++++++++++++++++++++++++"
stages = ["find", "ci", "vertex_coords", "rac_acts",
          "rac_inacts", "rho_acts",
          "rho_inacts", "tot_forces",
          "rgtp_forces", "edge_forces", "cyto_forces",
          "kdgtps_rac", "kgtps_rho", "kdgtps_rho", "complete"]

with open(test_file, 'r') as rf:
    lines = rf.readlines()

py_dat = [[], []]
step_dat = dict()
ci = None
stage = 0
for line in lines:
    line = line.strip()
    stage_string = stages[stage]
    if stage_string == "find" and delimiter in line:
        r = delimiter in line
        stage = (stage + 1) % len(stages)
        continue
    elif stage_string == "ci" and "ci" in line:
        r = stage_string in line
        dat_string = parse.parse("ci: {cell_index}", line).named["cell_index"]
        ci = eval(dat_string)
        step_dat[stage_string] = ci
        stage = (stage + 1) % len(stages)
        continue
    elif stage_string == "complete" and delimiter in line:
        r = delimiter in line
        py_dat[ci].append(copy.deepcopy(step_dat))
        step_dat = dict()
        ci = None
        stage = (stage + 1) % len(stages)
        continue
    elif stage_string in line:
        r = stage_string in line
        parse_string = "{}: {{{}}}".format(stage_string, stage_string)
        dat_string = eval(parse.parse(parse_string, line).named[stage_string])
        step_dat[stage_string] = np.array(dat_string)
        stage = (stage + 1) % len(stages)
        continue

world_history = None
snapshots = []
file_name = "../output/history_n_cells.cbor"
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

# deltas = dict()
# for label in labels:
#     print(label)
#     if label in ['vertex_coords', 'edge_forces', 'rgtpase_forces', 'cyto_forces', 'uivs']:
#         deltas[label] = np.linalg.norm(rust_dat[label] - py_dat[label], axis=2)
#     else:
#         deltas[label] = np.abs(rust_dat[label] - py_dat[label])
#


