import numpy as np
import parse
import matplotlib.pyplot as plt

data = []
last_was_data_line = False
with open("../output/console_output.txt", 'r') as f:
    for line in f.readlines():
        if line[:11] == "AsyncWriter":
            break
        if line[:13] == "ci_vi_oci_ovi":
            data.append(line)
            last_was_data_line = True
        elif last_was_data_line:
            data[-1] += line
            last_was_data_line = False

print("num lines found: {}".format(len(data)))
a0 = "ci_vi_oci_ovi:({ci},{vi},{oci},{ovi}),"
a = "zero_at:{zero_at},"
b = "adh_break:{adh_break},adh_rest:{adh_rest},"
c = "vc_mag:{vc_mag},adh_strain:{adh_strain},"
d = "smooth_factor:{smooth_factor}"
parse_string = a0 + a + b + c + d
print(parse_string)

keys = []
vc_mags = []
adh_strains = []
smooth_factors = []
zero_at = 0.0
adh_break = 0.0
adh_rest = 0.0
for ix, line in enumerate(data):
    line = line.strip()
    line = "".join(line.split())
    r = parse.parse(parse_string, line)
    vc_mags.append(float(r["vc_mag"]))
    adh_strains.append(float(r["adh_strain"]))
    zero_at = float(r["zero_at"])
    adh_break = float(r["adh_break"])
    adh_rest = float(r["adh_rest"])
    smooth_factors.append(float(r["smooth_factor"]))
    keys.append((int(r["ci"]), int(r["vi"]), int(r["oci"]), int(r["ovi"])))

set_keys = set(keys)
data = dict([(k, list()) for k in set_keys])
for k in set_keys:
    for i, x in enumerate(keys):
        if x == k:
            data[k].append([vc_mags[i], adh_strains[i], smooth_factors[i]])

for k in data.keys():
    data[k] = np.array(data[k])
print(data.keys())

d = data[(0, 0, 1, 8)]
changes = []
for x, y in zip(d[1:, 0], d[:-1, 0]):
    if np.abs(x - y) < 1e-3:
        changes.append(np.nan)
    elif x > y:
        changes.append(1)
    else:
        changes.append(-1)
plt.plot(np.arange(len(changes)), changes)
# # plt.plot(vc_mags, adh_strains, ls="", marker=".")
# begin = 0
# end = 10000
# plt.plot(vc_mags[begin:end], adh_strains[begin:end],
#          marker=".", ls="",
#          color="b", alpha=0.1)
# plt.plot(vc_mags[begin:end], smooth_factors[begin:end],
#          marker=".", ls="",
#          color="r", alpha=0.1)
