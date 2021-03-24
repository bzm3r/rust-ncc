import numpy as np
import parse
import matplotlib.pyplot as plt

data = []
last_was_data_line = False
with open("../output/console_output.txt", 'r') as f:
    for line in f.readlines():
        if line[:11] == "AsyncWriter":
            break
        if line[:7] == "zero_at":
            data.append(line)
            last_was_data_line = True
        elif last_was_data_line:
            data[-1] += line
            last_was_data_line = False

print("num lines found: {}".format(len(data)))
vc_mags = []
adh_strains = []
zero_at = 0.0
adh_max = 0.0
adh_delta_break = 0.0
adh_rest = 0.0
for ix, line in enumerate(data):
    line = line.replace("\n", "")
    r = parse.parse("zero_at: {zero_at}, adh_max: {adh_max}, adh_delta_break: {adh_delta_break}, adh_rest: {adh_rest}, vc_mag: {vc_mag}, adh_strain: {adh_strain}", line)
    vc_mags.append(float(r["vc_mag"]))
    adh_strains.append(float(r["adh_strain"]))
    zero_at = r["zero_at"]
    adh_max = r["adh_max"]
    adh_delta_break = r["adh_delta_break"]
    adh_rest = r["adh_rest"]

# plt.plot(vc_mags, adh_strains, ls="", marker=".")
begin = 0
end = 100000
plt.plot(vc_mags[begin:end], adh_strains[begin:end], marker=".", ls="")
# for ix in range(begin, end):
#     plt.annotate("{}".format(ix), (vc_mags[ix], adh_strains[ix]))
#
