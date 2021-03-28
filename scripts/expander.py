import os
import json
import copy
import numpy as np

exp_dir = "./experiments"
exp = "two_cell"

with open(os.path.join(exp_dir, exp + ".json")) as f:
    exp_dict = json.load(f)

cleanup = True
array_opt_labels = ["cil_mag", "coa_mag", "adh_mag", "cal_mag", "one_at", "zero_at", "too_close_dist"]
array_opts = dict()
expected_num_combinations = np.prod

for label in array_opt_labels:
    if label in exp_dict.keys():
        data = exp_dict[label]
        if type(data) == list:
            array_opts[label] = data
        else:
            array_opts[label] = [data]

found_labels = [label for label in array_opt_labels if label in
                array_opts.keys()]
total_combos = np.prod([len(array_opts[key]) for key in array_opts.keys()])

combinations = []
for label in found_labels:
    if len(combinations) == 0:
        combinations = [list(array_opts[label]) for label in array_opts.keys()]



# combinations = []
# for label in found_labels:
#     init_combinations = copy.deepcopy(combinations)
#     new_combinations = []
#     for v in array_opts[label]:
#         if len(init_combinations) == 0:
#             new_combinations.append([(label, v)])
#         else:
#             for c in init_combinations:
#                 new_combinations.append(copy.deepcopy(init_combinations) + [(label, v)])
#     combinations += copy.deepcopy(new_combinations)
#
