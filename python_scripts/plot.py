import matplotlib.pyplot as plt
import json
import numpy as np
import json
import cbor2


output = None

file_name = "history_cal_test.cbor"

with open(file_name, mode='rb') as sf:
    output = cbor2.load(sf)

tsteps = [o[0] for o in output]
frequency = tsteps[1] - tsteps[0]
state_recs = [o[1] for o in output]
interactions = [rec["interactions"] for rec in state_recs]
x_cals_0_4 = [inters[0]["x_cals"][4] for inters in interactions]
x_cals_1_12 = [inters[1]["x_cals"][12] for inters in interactions]


plt.plot(tsteps, x_cals_0_4, color="black", marker=".")
plt.plot(tsteps, x_cals_1_12, color="green", marker=".")
