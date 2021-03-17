import numpy as np
import matplotlib.pyplot as plt
import os
import cbor2
from matplotlib import animation
from paint_opts import *
import get_comp as cd
import get_cbor as cb
import copy

class SimulationData:
    out_dir = None
    file_name = None
    cbor_file_path = None
    mp4_file_path = None
    dat_file_path = None
    tpoints = None

    world_info = None
    header = None

    poly_per_c_per_s = None
    centroids_per_c_per_s = None
    uivs_per_c_per_s = None
    uovs_per_c_per_s = None

    rac_acts_per_c_per_s = None
    rac_inacts_per_c_per_s = None
    rac_act_arrows_per_c_per_s = None

    rho_acts_per_c_per_s = None
    rho_inacts_per_c_per_s = None
    rho_act_arrows_per_c_per_s = None

    kgtps_rac_per_c_per_s = None
    kdgtps_rac_per_c_per_s = None
    kgtps_rho_per_c_per_s = None
    kdgtps_rho_per_c_per_s = None
    rac_act_net_fluxes_per_c_per_s = None
    rac_inact_net_fluxes_per_c_per_s = None
    rho_act_net_fluxes_per_c_per_s = None
    rho_inact_net_fluxes_per_c_per_s = None
    x_tens_per_c_per_s = None

    x_cils_per_c_per_s = None
    x_coas_per_c_per_s = None
    x_adhs_per_c_per_s = None

    edge_strains_per_c_per_s = None
    rgtp_forces_per_c_per_s = None
    edge_forces_per_c_per_s = None
    cyto_forces_per_c_per_s = None
    sum_forces_per_c_per_s = None
    avg_tens_strain_per_c_per_s = None

    snap_ix = None
    default_xlim = None
    default_ylim = None
    default_bbox_lim = None
    ani_opts = None
    mp4_file_name_header = None
    fig_probe = None
    ax_probe = None
    fig_ani = None
    ax_ani = None
    snap_period = None
    char_t = None

    def generate_header_from_world_info(self):
        all_params = dict()
        all_params.update(self.world_info["char_quants"])
        world_params = self.world_info["world_params"]
        vertex_eta = world_params["vertex_eta"]
        all_params["vertex_eta"] = vertex_eta
        inter_params = world_params["interactions"]
        phys_params = inter_params["phys_contact"]
        range_params = phys_params["range"]
        all_params.update(range_params)
        all_params["cil_mag"] = phys_params["cil_mag"]
        coa_params = inter_params["coa"]
        for key in coa_params.keys():
            all_params["coa_" + key] = coa_params[key]
        all_params.update(self.world_info["cell_params"][0])
        self.header = copy.deepcopy(all_params)

    def load_rust_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.cbor_file_path = self.file_name + ".cbor"
        self.mp4_file_name_header = self.file_name + "_M=r"

        cbor_files = \
            [f for f in os.listdir(self.out_dir)
             if os.path.isfile(os.path.join(self.out_dir, f))
             and os.path.splitext(f)[1] == ".cbor"]

        found_wanted = False
        for fn in cbor_files:
            if self.cbor_file_path == fn:
                found_wanted = True
                break

        if not found_wanted:
            raise Exception(
                "Error: could not find requested file {} in dir {} with "
                "contents: {}".format(self.cbor_file_path, out_dir, cbor_files))
        self.cbor_file_path = os.path.join(self.out_dir, self.cbor_file_path)

        snapshots = []
        with open(self.cbor_file_path, mode='rb') as sf:
            world_info = cbor2.load(sf)
            while True:
                try:
                    snapshots += cbor2.load(sf)
                except EOFError:
                    break

        print("load_rust_dat | file_name: {} | snapshots found: {}"
              .format(file_name, len(snapshots)))
        self.world_info = world_info
        self.generate_header_from_world_info()
        self.char_t = world_info["char_quants"]["t"]
        self.tpoints = [s["cells"][0]["tpoint"] * self.char_t for s in
                        snapshots]
        data = [s["cells"] for s in snapshots]
        self.snap_period = world_info["snap_period"]

        self.poly_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'poly'], data)
        self.centroids_per_c_per_s = np.array(
            [[np.average(poly, axis=0) for poly in poly_per_c] for
             poly_per_c in
             self.poly_per_c_per_s])
        self.uivs_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'geom', 'unit_in_vecs'],
                                      data)
        self.uovs_per_c_per_s = -1 * self.uivs_per_c_per_s
        self.rac_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rac_acts'], data)
        self.rac_inacts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rac_inacts'], data)
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_acts'], data)
        self.rho_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_inacts'], data)
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s

        self.x_cils_per_c_per_s = \
            cb.extract_scalars_from_data(['interactions', 'x_cils'], data)
        self.x_coas_per_c_per_s = \
            cb.extract_scalars_from_data(['interactions', 'x_coas'], data)
        self.x_adhs_per_c_per_s = \
            cb.extract_p2ds_from_data(['interactions', 'x_adhs'], data)

        self.kgtps_rac_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'kgtps_rac'], data)
        self.kdgtps_rac_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'kdgtps_rac'], data)
        self.kgtps_rho_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'kgtps_rho'], data)
        self.kdgtps_rho_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'kdgtps_rho'], data)

        self.rac_act_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'rac_act_net_fluxes'],
                                         data)
        self.rac_inact_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['chem',
                                          'rac_inact_net_fluxes'], data)
        self.rho_act_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'rho_act_net_fluxes'],
                                         data)
        self.rho_inact_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['chem',
                                          'rho_inact_net_fluxes'], data)

        self.x_tens_per_c_per_s = \
            cb.extract_scalars_from_data(['chem', 'x_tens'], data)
        self.edge_strains_per_c_per_s = \
            cb.extract_scalars_from_data(['mech', 'edge_strains'], data)
        self.rgtp_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['mech', 'rgtp_forces'], data)
        self.edge_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['mech', 'edge_forces'], data)
        self.cyto_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['mech', 'cyto_forces'], data)
        self.sum_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['mech', 'sum_forces'], data)
        self.avg_tens_strain_per_c_per_s = \
            cb.extract_scalars_from_data(['mech', 'avg_tens_strain'],
                                         data)

    def load_py_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.dat_file_path = self.file_name + ".dat"
        self.mp4_file_name_header = self.file_name + "_M=p"

        raw_out = cd.read_save_file(self.out_dir, self.dat_file_path)
        data_per_c_per_s = cd.get_data_per_c_per_s(raw_out)
        print("load_py_dat | file_name: {} | snapshots found: {}"
              .format(file_name, len(data_per_c_per_s)))
        self.header = raw_out["header"]
        num_int_steps = raw_out["header"]["num_int_steps"]
        self.char_t = raw_out["header"]["t"]
        self.tpoints = [cells[0]["tpoint"] * self.char_t for cells in
                        data_per_c_per_s]
        self.snap_period = num_int_steps

        self.poly_per_c_per_s = \
            np.array([[cell["poly"] for cell in cells]
                      for cells in data_per_c_per_s])
        self.centroids_per_c_per_s = np.array(
            [[np.average(poly, axis=0) for poly in poly_per_c] for
             poly_per_c in
             self.poly_per_c_per_s])
        self.uivs_per_c_per_s = \
            np.array([[snap["uivs"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.uovs_per_c_per_s = -1 * self.uivs_per_c_per_s
        self.rac_acts_per_c_per_s = \
            np.array([[snap["rac_acts"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rac_inacts_per_c_per_s = \
            np.array([[snap["rac_inacts"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = np.array([[snap["rho_acts"] for snap in
                                               snaps_per_c]
                                              for snaps_per_c in
                                              data_per_c_per_s])
        self.rho_inacts_per_c_per_s = \
            np.array([[snap["rho_acts"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s

        self.x_cils_per_c_per_s = \
            np.array([[snap["x_cils"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.x_coas_per_c_per_s = \
            np.array([[snap["x_coas"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])

        self.kgtps_rac_per_c_per_s = \
            np.array([[snap["kgtps_rac"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.kdgtps_rac_per_c_per_s = \
            np.array([[snap["kdgtps_rac"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.kgtps_rho_per_c_per_s = \
            np.array([[snap["kgtps_rho"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.kdgtps_rho_per_c_per_s = \
            np.array([[snap["kdgtps_rho"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])

        self.rac_act_net_fluxes_per_c_per_s = \
            self.rac_inact_net_fluxes_per_c_per_s = \
            np.array([[snap["rac_act_net_fluxes"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rac_inact_net_fluxes_per_c_per_s = \
            np.array([[snap["rac_inact_net_fluxes"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])

        self.rho_act_net_fluxes_per_c_per_s = \
            self.rho_inact_net_fluxes_per_c_per_s = \
            np.array([[snap["rho_act_net_fluxes"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rho_inact_net_fluxes_per_c_per_s = \
            np.array([[snap["rho_inact_net_fluxes"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])

        self.x_tens_per_c_per_s = \
            np.array([[snap["x_tens"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.edge_strains_per_c_per_s = \
            np.array([[snap["edge_strains"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rgtp_forces_per_c_per_s = \
            np.array([[snap["rgtp_forces"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.edge_forces_per_c_per_s = \
            np.array([[snap["edge_forces"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.cyto_forces_per_c_per_s = \
            np.array([[snap["cyto_forces"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.sum_forces_per_c_per_s = \
            np.array([[snap["sum_forces"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.avg_tens_strain_per_c_per_s = \
            np.array([[snap["avg_tens_strain"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])