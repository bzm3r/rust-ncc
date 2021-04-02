import numpy as np
import matplotlib.pyplot as plt


# ====================================
# function to calculate centroid info
def calculate_centroid_info(centroids_per_c_per_s):
    # deltas_per_c_per_s should have shape:
    # [num_snapshots - 1, num_cells, 2]
    delta_per_c_per_s = centroids_per_c_per_s[1:] - centroids_per_c_per_s[:-1]

    # calculate the distance per snapshot, per cell, will have shape:
    # [num_snapshots - 1, num_cells]
    dist_per_c_per_s = np.linalg.norm(delta_per_c_per_s, axis=2)

    # ================================

    # avg_centroid_per_s should have shape [num_snapshots, 2]
    avg_centroid_per_s = np.average(centroids_per_c_per_s, axis=1)

    # delta_group_centroid_per_s will have shape: [num_snapshots - 1, 2]
    delta_group_centroid_per_s = avg_centroid_per_s[1:] - avg_centroid_per_s[
                                                          :-1]

    # calculate the distance per snapshot, will have shape: [num_snapshots - 1]
    dist_group_centroid_per_s = np.linalg.norm(delta_group_centroid_per_s,
                                               axis=1)

    return dist_per_c_per_s, dist_group_centroid_per_s


# ===============================
# an example assuming only one experiment, one seed, with random info
num_cells = 16
num_snapshots = 10

# create some random data
centroids_per_c_per_s = np.random.rand(num_snapshots, num_cells, 2)

# first get calculate centroid info
dist_per_c_per_s, dist_group_centroid_per_s = \
    calculate_centroid_info(centroids_per_c_per_s)

# then do plots
path_dist_cell_0 = np.sum(dist_per_c_per_s[:, 0])
print("path_dist_cell_0: {}".format(path_dist_cell_0))
end_to_end_dist_cell_0 = np.linalg.norm(centroids_per_c_per_s[-1, 0] -
                                        centroids_per_c_per_s[0, 0])
print("end_to_end_dist_cell_0: {}".format(end_to_end_dist_cell_0))
plt.plot(centroids_per_c_per_s[:, 0, 0], centroids_per_c_per_s[:, 0, 1],
         color="b")

path_dist_cell_0 = np.sum(dist_group_centroid_per_s[:])
print("path_dist_group: {}".format(path_dist_cell_0))
end_to_end_dist_cell_0 = np.linalg.norm(avg_centroid_per_s[-1] -
                                        avg_centroid_per_s[0])
print("end_to_end_dist_group: {}".format(end_to_end_dist_cell_0))
plt.plot(avg_centroid_per_s[:, 0], avg_centroid_per_s[:, 1], color="r")

# persistence = end_to_end / path_dist

# =============================================================
# how to get the data from the experiments?
root_dir = os.getcwd()
out_dir = os.path.join(root_dir, "output")
exp_jsons = ["adh_0, adh_5"]

seeded_experiments_per_json = []
# maybe exp_jsons is ["adh_0", "adh_5"]
for exp_json in exp_json:
    exp_path = os.path.join(root_dir, "experiments", "{}.json".format(exp_json))
    with open(exp_path) as f:
        json_str = f.read()
    exp_dict = orjson.loads(json_str)
    # get seeds, and file name per seed, for experiment defined in exp_json
    seeds, file_names = determine_file_names(exp_json, exp_dict)

    seeded_experiments_this_json = []
    # maybe seeds are [7, 77, 77]
    for seed, file_name in zip(seeds, file_names):
        rust_dat = SimulationData()
        rust_dat.load_rust_dat(out_dir, file_name)
        rust_dat.seed = seed
        seeded_experiments_this_json.append(rust_dat)
    seeded_experiments_per_json.append(seeded_experiments_this_json)

centroids_per_c_per_s_per_seed_per_experiment = \
    [[rust_dat.centroids_per_c_per_s
      for rust_dat in rust_dats_per_seed] for rust_dats_per_seed \
     in rust_dats_per_seed_per_json]

# now we'd have centroid info for each seed, in each experiment
dist_group_centroid_per_s_per_seed_per_experiment = \
    [[calculate_centroid_info(centroids_per_c_per_s) for
      centroids_per_c_per_s in
      centroids_per_seed] for
     centroids_per_seed in
     centroids_per_c_per_s_per_seed_per_experiment]

# =========================================================================
