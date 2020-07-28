import os
import numpy as np


labels = ["vertex_coords", "rac_acts",
          "rac_inacts", "rho_acts",
          "rho_inacts", "uivs",
          "rgtp_forces", "edge_strains", "edge_forces",
          "avg_tens_strain", "cyto_forces"]

def process_lines(fp):
    with open(fp, 'r') as rf:
        lines = rf.readlines()

    results = dict([(label, list()) for label in labels])
    for line in lines:
        for label in labels:
            n = len(label)
            if line[:n] == label:
                dat_str = line[n+2:].strip()
                results[label].append(np.array(eval(dat_str)))
    
    for label in labels:
        results[label] = np.array(results[label])
        
    return results

test_file = "out_euler.txt"

rust_dir = "C:\\Users\\bhmer\\Desktop\\rust-ncc\\"
py_dir = "C:\\Users\\bhmer\\Desktop\\numba-ncc\\"

rust_dat = process_lines(os.path.join(rust_dir, test_file))
py_dat = process_lines(os.path.join(py_dir, test_file))

deltas = dict()
for label in labels:
    print(label)
    if label in ['vertex_coords', 'edge_forces', 'rgtpase_forces', 'cyto_forces', 'uivs']:
        deltas[label] = np.linalg.norm(rust_dat[label] - py_dat[label], axis=2)
    else:
        deltas[label] = np.abs(rust_dat[label] - py_dat[label])



