import numpy as np
import matplotlib.pyplot as plt
import os
import cbor2
from matplotlib import animation
from paint_opts import *
import get_comp as cd
import get_cbor as cb
from utils import *
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
    tag = None

    poly_per_c_per_s = None
    centroids_per_c_per_s = None
    uivs_per_c_per_s = None
    uovs_per_c_per_s = None

    rac_acts_per_c_per_s = None
    rac_inacts_per_c_per_s = None
    rac_acts_arrows_per_c_per_s = None

    rho_acts_per_c_per_s = None
    rho_inacts_per_c_per_s = None
    rho_acts_arrows_per_c_per_s = None
    rho_inacts_arrows_per_c_per_s = None

    rho_inacts_arrow_group = None
    rho_acts_arrow_group = None
    rgtps_arrow_group = None

    kgtps_rac_per_c_per_s = None
    kgtps_rac_arrows_per_c_per_s = None
    kgtps_rac_arrow_group = None

    kdgtps_rac_per_c_per_s = None
    kdgtps_rac_arrows_per_c_per_s = None
    kdgtps_rac_arrow_group = None

    kgtps_rho_per_c_per_s = None
    kgtps_rho_arrows_per_c_per_s = None
    kgtps_rho_arrow_group = None

    kdgtps_rho_per_c_per_s = None
    kdgtps_rho_arrows_per_c_per_s = None
    kdgtps_rho_arrow_group = None

    rac_act_net_fluxes_per_c_per_s = None
    rac_inact_net_fluxes_per_c_per_s = None
    rho_act_net_fluxes_per_c_per_s = None
    rho_act_net_fluxes_arrows_per_c_per_s = None
    rho_act_net_fluxes_arrow_group = None
    rho_inact_net_fluxes_per_c_per_s = None
    x_tens_per_c_per_s = None
    x_tens_arrows_per_c_per_s = None
    x_tens_arrow_group = None

    x_cils_per_c_per_s = None
    x_cals_per_c_per_s = None
    x_cils_arrows_per_c_per_s = None
    x_cals_arrows_per_c_per_s = None
    x_cils_arrow_group = None
    x_cals_arrow_group = None
    x_coas_per_c_per_s = None
    x_coas_arrows_per_c_per_s = None
    x_coas_arrow_group = None
    x_adhs_per_c_per_s = None

    edge_strains_per_c_per_s = None
    rgtp_forces_per_c_per_s = None
    rgtp_forces_arrow_group = None
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
        if coa_params != None:
            for key in coa_params.keys():
                all_params["coa_" + key] = coa_params[key]
        all_params.update(self.world_info["cell_params"][0])

        all_params["init_rac"] = all_params["init_rac"]["active"]
        all_params["init_rho"] = all_params["init_rho"]["active"]
        self.header = copy.deepcopy(all_params)

    def load_animation_arrows(self):
        self.uovs_per_c_per_s = -1 * self.uivs_per_c_per_s

        self.rac_acts_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.rho_inacts_arrows_per_c_per_s = \
            self.rho_inacts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.rgtps_arrow_group = [(100.0, "b",
                                   self.rac_acts_arrows_per_c_per_s),
                                  (
                                      100.0, "r",
                                      self.rho_acts_arrows_per_c_per_s)]
        self.rho_acts_arrow_group = [(500.0, "r",
                                      self.rho_acts_arrows_per_c_per_s)]
        self.rho_inacts_arrow_group = [(100.0, "r",
                                        self.rho_inacts_arrows_per_c_per_s)]

        self.x_cils_arrows_per_c_per_s = \
            self.x_cils_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.x_cils_arrow_group = [(0.5, "r", self.x_cils_arrows_per_c_per_s)]

        if self.tag == "rust":
            self.x_cals_arrows_per_c_per_s = \
                self.x_cals_per_c_per_s[:, :, :, np.newaxis] * \
                self.uovs_per_c_per_s
            self.x_cals_arrow_group = [
                (0.1, "r", self.x_cals_arrows_per_c_per_s)]

        self.x_coas_arrows_per_c_per_s = \
            self.x_coas_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.x_coas_arrow_group = [(3.0, "b", self.x_coas_arrows_per_c_per_s)]

        self.kgtps_rho_arrows_per_c_per_s = \
            self.kgtps_rho_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.kgtps_rho_arrow_group = [(100.0, "r",
                                       self.kgtps_rho_arrows_per_c_per_s)]

        self.kdgtps_rho_arrows_per_c_per_s = \
            self.kdgtps_rho_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.kdgtps_rho_arrow_group = [(100.0, "r",
                                        self.kdgtps_rho_arrows_per_c_per_s)]

        self.kgtps_rac_arrows_per_c_per_s = \
            self.kgtps_rac_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.kgtps_rac_arrow_group = [(100.0, "b",
                                       self.kgtps_rac_arrows_per_c_per_s)]

        self.kdgtps_rac_arrows_per_c_per_s = \
            self.kdgtps_rac_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.kdgtps_rac_arrow_group = [(100.0, "b",
                                        self.kdgtps_rac_arrows_per_c_per_s)]

        self.rho_act_net_fluxes_arrows_per_c_per_s = \
            self.rho_act_net_fluxes_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.rho_act_net_fluxes_arrow_group = [(100.0, "b",
                                                self.rho_act_net_fluxes_arrows_per_c_per_s)]

        is_rac_force = np.zeros_like(self.rgtp_forces_per_c_per_s)
        is_rho_force = np.zeros_like(self.rgtp_forces_per_c_per_s)

        ones2d = np.ones(2)
        for s_ix in range(is_rac_force.shape[0]):
            for c_ix in range(is_rac_force.shape[1]):
                for v_ix in range(is_rac_force.shape[2]):
                    uov = self.uovs_per_c_per_s[s_ix, c_ix, v_ix]
                    rgtp_f = self.rgtp_forces_per_c_per_s[s_ix, c_ix, v_ix]
                    dp = np.dot(uov, rgtp_f)

                    if abs(dp) > 1e-4:
                        if dp < 0.0:
                            is_rho_force[s_ix, c_ix, v_ix, :] = ones2d
                        else:
                            is_rac_force[s_ix, c_ix, v_ix, :] = ones2d

        self.rgtp_forces_arrow_group = [(0.1, "orange",
                                         self.rgtp_forces_per_c_per_s * is_rho_force),
                                        (0.1, "green",
                                         self.rgtp_forces_per_c_per_s * is_rac_force)]

    def load_rust_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.cbor_file_path = self.file_name + ".cbor"
        self.mp4_file_name_header = self.file_name + "_M=r"
        self.tag = "rust"

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

        self.rac_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rac_acts'], data)
        self.rac_inacts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rac_inacts'], data)
        self.rho_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_acts'], data)
        self.rho_inacts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_inacts'], data)

        self.x_cils_per_c_per_s = \
            cb.extract_scalars_from_data(['interactions', 'x_cils'], data)
        self.x_cals_per_c_per_s = \
            cb.extract_scalars_from_data(['interactions', 'x_cals'], data)
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

        self.load_animation_arrows()

    def load_py_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.dat_file_path = self.file_name + ".dat"
        self.mp4_file_name_header = self.file_name + "_M=p"
        self.tag = "py"

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
        self.rho_acts_per_c_per_s = np.array([[snap["rho_acts"] for snap in
                                               snaps_per_c]
                                              for snaps_per_c in
                                              data_per_c_per_s])
        self.rho_inacts_per_c_per_s = \
            np.array([[snap["rho_inacts"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])

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

        self.load_animation_arrows()

    def animate(self, vec_ani_opts, ty):
        self.default_xlim = [-40, 200]
        self.default_ylim = [-40, 200]
        self.default_bbox_lim = \
            [self.default_xlim[1] - self.default_xlim[0],
             self.default_ylim[1] - self.default_ylim[0]]
        self.fig_ani, self.ax_ani = plt.subplots()
        self.ax_ani.set_aspect('equal')
        self.ax_ani.set_xlim(self.default_xlim)
        self.ax_ani.set_ylim(self.default_ylim)
        sanitized_tpoints = sanitize_tpoints(self.tpoints)
        frame_ixs = [x[0] for x in sanitized_tpoints]
        for ani_opts in vec_ani_opts:
            self.snap_ix = 0
            self.ani_opts = ani_opts
            # Set up formatting for the movie files
            writer = animation.writers['ffmpeg'](fps=10,
                                                 metadata=dict(artist='Me'),
                                                 bitrate=1800)
            cell_ani = animation.FuncAnimation(self.fig_ani,
                                               self.paint_animation,
                                               frames=frame_ixs,
                                               fargs=(self.ax_ani, ty),
                                               interval=1, blit=True)
            if len(ty) != 0:
                ty_tag = "_{}".format(ty)
            else:
                ty_tag = ""

            ani_file_path = self.mp4_file_name_header + ani_opts.description(

            ) + ty_tag + ".mp4"
            ani_save_path = os.path.join(self.out_dir, ani_file_path)
            cell_ani.save(ani_save_path, writer=writer)

    def paint_animation(self, snap_ix, ax, ty):
        ax.cla()
        ax.set_aspect("equal")
        print("painting snapshot: {}".format(snap_ix))
        if not self.ani_opts.follow_group:
            g_center = np.average(self.centroids_per_c_per_s[snap_ix],
                                  axis=0)
            (xmin, xmax) = [g_center[0] - self.default_bbox_lim[0] * 0.5,
                            g_center[0] + self.default_bbox_lim[0] * 0.5]
            (ymin, ymax) = [g_center[1] - self.default_bbox_lim[1] * 0.5,
                            g_center[1] + self.default_bbox_lim[1] * 0.5]
            bbox = \
                np.array([
                    [xmin, ymin],
                    [xmin, ymax],
                    [xmax, ymax],
                    [xmax, ymin],
                    [xmin, ymin]
                ])
            ax.plot(bbox[:, 0], bbox[:, 1], color=(0.0, 0.0, 0.0, 0.0))

        self.paint_cells(snap_ix, ax, ty)
        ax.relim()

        ax.set_title(
            "tstep {} (snapshot {})".format(
                self.tpoints[snap_ix], snap_ix
            )
        )
        return ax.get_children()

    def paint_cells(self, snap_ix, ax, ty):
        pls = self.ani_opts.poly_line_style
        for (ci, poly) in enumerate(self.poly_per_c_per_s[snap_ix]):
            for vix in range(16):
                ax.plot([poly[vix, 0], poly[(vix + 1) % 16, 0]],
                        [poly[vix, 1], poly[(vix + 1) % 16, 1]],
                        color="k", marker=".", markersize="0.5",
                        linestyle=pls)
                if self.ani_opts.label_verts:
                    ax.annotate(str(vix), (poly[vix, 0], poly[vix, 1]))

            c_centers = self.centroids_per_c_per_s[:snap_ix, ci]
            if self.ani_opts.show_trails:
                ax.plot(c_centers[:, 0], c_centers[:, 1])
            if self.ani_opts.label_cells and snap_ix > 0:
                ax.annotate(str(ci), (c_centers[-1, 0],
                                      c_centers[-1, 1]))

        arrow_group = eval("self.{}_arrow_group".format(ty))

        for (scale, color, arrows_per_c_per_s) in arrow_group:
            for poly, arrows in zip(
                    self.poly_per_c_per_s[snap_ix],
                    arrows_per_c_per_s[snap_ix]
            ):
                for p, arrow in zip(poly, arrows):
                    ax.arrow(p[0], p[1],
                             scale * self.ani_opts.arrow_scale * arrow[0],
                             scale * self.ani_opts.arrow_scale * arrow[1],
                             color=color, linestyle=pls,
                             length_includes_head=True, head_width=0.0)
