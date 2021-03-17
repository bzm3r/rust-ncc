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

    def load_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.cbor_file_path = self.file_name + ".cbor"
        self.mp4_file_name_header = self.file_name + ".mp4"

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
            success = True
            while success:
                try:
                    snapshots += cbor2.load(sf)
                finally:
                    success = False

        print("cbor file: {} | snapshots found: {}"
              .format(file_name, len(snapshots)))
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
            cb.extract_scalars_from_data(['core', 'chem', 'kgtps_rac'], data)
        self.kdgtps_rac_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'kdgtps_rac'], data)
        self.kgtps_rho_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'kgtps_rho'], data)
        self.kdgtps_rho_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'kdgtps_rho'], data)

        self.rac_act_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'rac_act_net_fluxes'],
                                      data)
        self.rac_inact_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem',
                                       'rac_inact_net_fluxes'], data)
        self.rho_act_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'rho_act_net_fluxes'],
                                      data)
        self.rho_inact_net_fluxes_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem',
                                          'rho_inact_net_fluxes'], data)

        self.x_tens_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'chem', 'x_tens'], data)
        self.edge_strains_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'mech', 'x_tens'], data)
        self.rgtp_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'mech', 'x_tens'], data)
        self.edge_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'mech', 'x_tens'], data)
        self.cyto_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'mech', 'x_tens'], data)
        self.sum_forces_per_c_per_s = \
            cb.extract_p2ds_from_data(['core', 'mech', 'x_tens'], data)
        self.avg_tens_strain_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'mech', 'x_tens'], data)

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
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_acts'], data)
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.x_adhs_per_c_per_s = \
            cb.extract_p2ds_from_data(['interactions', 'x_adhs'], data)

    def load_py_dat(self, out_dir, file_name):
        self.out_dir = out_dir
        self.file_name = file_name
        self.dat_file_path = self.file_name + ".dat"
        self.mp4_file_name_header = self.file_name + "_M=p"

        raw_out = cd.read_save_file(self.out_dir, self.dat_file_path)
        data_per_c_per_s = cd.get_data_per_c_per_s(raw_out)
        print("load_py_dat | file_name: {} | snapshots found: {}"
              .format(file_name, len(data_per_c_per_s)))
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
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = \
            np.array([[snap["rho_acts"] for snap in snaps_per_c]
                      for snaps_per_c in data_per_c_per_s])
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.x_adhs_per_c_per_s = np.zeros_like(self.poly_per_c_per_s)

    def probe(self, ani_opts):
        self.ani_opts = ani_opts
        self.snap_ix = 0
        self.default_xlim = [-100, 200]
        self.default_ylim = [-100, 200]
        self.fig_probe, self.ax_probe = plt.subplots()
        self.ax_probe.set_aspect('equal')
        self.ax_probe.set_xlim(self.default_xlim)
        self.ax_probe.set_ylim(self.default_ylim)
        self.fig_probe.canvas.mpl_connect(
            'key_press_event',
            lambda event: self.on_press_probe(event)
        )
        plt.show()

    def animate(self, vec_ani_opts):
        self.default_xlim = [-40, 200]
        self.default_ylim = [-40, 200]
        self.default_bbox_lim = \
            [self.default_xlim[1] - self.default_xlim[0],
             self.default_ylim[1] - self.default_ylim[0]]
        self.fig_ani, self.ax_ani = plt.subplots()
        self.ax_ani.set_aspect('equal')
        self.ax_ani.set_xlim(self.default_xlim)
        self.ax_ani.set_ylim(self.default_ylim)
        frame_ixs = [n for n in range(int(len(self.tpoints)))]
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
                                               fargs=(self.ax_ani,),
                                               interval=1, blit=True)
            ani_file_path = self.mp4_file_name_header + ani_opts.description() + ".mp4"
            ani_save_path = os.path.join(self.out_dir, ani_file_path)
            cell_ani.save(ani_save_path, writer=writer)

    def paint_cells(self, snap_ix, ax):
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
            if self.ani_opts.label_cells:
                ax.annotate(str(ci), (c_centers[-1, 0], c_centers[-1, 1]))

        for poly, rac_act_arrows in zip(
                self.poly_per_c_per_s[snap_ix],
                self.rac_act_arrows_per_c_per_s[snap_ix]
        ):
            for p, rac_arrow in zip(poly, rac_act_arrows):
                ax.arrow(p[0], p[1], self.ani_opts.rgtp_scale * rac_arrow[0],
                         self.ani_opts.rgtp_scale * rac_arrow[1],
                         color="b",
                         length_includes_head=True, head_width=0.0)

        for poly, rho_act_arrows in zip(self.poly_per_c_per_s[snap_ix],
                                        self.rho_act_arrows_per_c_per_s[
                                            snap_ix]):
            for p, rho_arrow in zip(poly, rho_act_arrows):
                ax.arrow(p[0], p[1], self.ani_opts.rgtp_scale * rho_arrow[0],
                         self.ani_opts.rgtp_scale * rho_arrow[1],
                         color="r",
                         length_includes_head=True, head_width=0.0)

        for poly_ix, poly, adhs in zip(
                np.arange(0, len(self.poly_per_c_per_s[0])),
                self.poly_per_c_per_s[snap_ix],
                self.x_adhs_per_c_per_s[snap_ix]):
            if poly_ix == 0:
                adh_arrow_color = "magenta"
            else:
                adh_arrow_color = "cyan"
            for p, adh in zip(poly, adhs):
                ax.arrow(p[0], p[1], adh[0], adh[1], color=adh_arrow_color,
                         length_includes_head=True, head_width=1.0)

    def paint_probe(self, delta):
        old_xlim = self.ax_probe.get_xlim()
        old_ylim = self.ax_probe.get_ylim()
        self.ax_probe.cla()
        self.ax_probe.set_xlim(old_xlim)
        self.ax_probe.set_ylim(old_ylim)

        self.paint_cells(self.snap_ix, self.ax_probe)

        self.ax_probe.set_title(
            "tstep {} (snapshot {})".format(
                self.tpoints[self.snap_ix], self.snap_ix
            )
        )
        self.snap_ix = (self.snap_ix + delta) % len(self.tpoints)
        plt.show()

    def on_press_probe(self, event):
        if event.key == 'x':
            self.paint_probe(1)
        elif event.key == 'z':
            self.paint_probe(-1)
        if event.key == 'c':
            self.paint_probe(-5)
        elif event.key == 'v':
            self.paint_probe(5)
        elif event.key == 'n':
            self.paint_probe(-10)
        elif event.key == 'm':
            self.paint_probe(10)
        elif event.key == 'r':
            self.ax_probe.set_aspect('equal')
            self.ax_probe.set_xlim(self.default_xlim)
            self.ax_probe.set_ylim(self.default_ylim)
        self.fig_probe.canvas.draw()

    def paint_animation(self, snap_ix, ax):
        ax.cla()
        ax.set_aspect("equal")

        if self.ani_opts.follow_group:
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

        self.paint_cells(snap_ix, ax)
        ax.relim()

        ax.set_title(
            "tstep {} (snapshot {})".format(
                self.tpoints[snap_ix], snap_ix
            )
        )
        return ax.get_children()


def find_common_ts(ixs_ts_per_sim):
    if len(ixs_ts_per_sim) == 0:
        return []
    elif len(ixs_ts_per_sim) == 1:
        return ([t for (ix, t) in ixs_ts_per_sim[0]],
                [[ix for (ix, t) in ixs_ts] for ixs_ts in ixs_ts_per_sim])
    else:
        common = []
        ixs_per_list = [list() for _ in range(len(ixs_ts_per_sim))]
        num_commons = 0
        ixs_ts_per_sim = sorted(ixs_ts_per_sim, key=lambda x: len(x))
        xs = ixs_ts_per_sim[0]
        num_xs = len(xs)
        min_ixs = [0 for _ in ixs_ts_per_sim[1:]]
        for (x_snap_ix, x) in xs:
            ix_per_list = [x_snap_ix]
            if num_commons < num_xs:
                inside_all = True
                for (list_ix, ys) in enumerate(ixs_ts_per_sim[1:]):
                    if inside_all:
                        check = False
                        for (ix_y, (y_snap_ix, y)) in enumerate(ys):
                            if abs((x - y)) < 1e-4:
                                check = True
                                ix_per_list.append(y_snap_ix)
                                break
                        inside_all = check and inside_all
                    else:
                        break
                if inside_all:
                    common.append(x)
                    if len(ix_per_list) != len(ixs_per_list):
                        raise Exception(
                            "len(ix_per_list) = {} != {} = len(ixs_per_list)"
                            .format(len(ix_per_list), len(ixs_per_list)))
                    for (ix, ixs) in zip(ix_per_list, ixs_per_list):
                        ixs.append(ix)

                    num_commons += 1

        return common, ixs_per_list


def sanitize_tpoints(tpoints):
    tpoints = np.array(tpoints)
    deltas = tpoints[1:] - tpoints[:-1]
    max_delta = np.max(deltas)
    is_max = np.flatnonzero(np.abs(deltas - max_delta) < 1e-4) + 1
    sanitized_tpoints = [tpoints[0]] + tpoints[is_max].tolist()
    snap_ixs = [0] + is_max.tolist()
    return list(zip(snap_ixs, sanitized_tpoints))


class SharedSimData:
    def __init__(self, out_dir, sim_dats, poly_line_styles, mp4_file_name):
        self.out_dir = out_dir
        self.sim_dats = sim_dats
        self.mp4_file_name = mp4_file_name
        self.poly_line_styles = poly_line_styles

        self.vert_plot_ix = 0
        self.curr_inner_ix = 0
        self.plot_x_max = 0
        self.dat_group_ix = 0
        self.cell_plot_ix = 0
        self.num_cells = 0
        self.dat_groups = 0
        self.active_dat_group_ix = 0
        self.num_label_groups = 0

        self.common_ts, self.snap_ixs_per_sim = self.get_common_snaps()
        self.snap_ix = 0
        self.default_xlim = [-40, 200]
        self.default_ylim = [-40, 200]
        self.default_bbox_lim = \
            [self.default_xlim[1] - self.default_xlim[0],
             self.default_ylim[1] - self.default_ylim[0]]
        self.fig, self.ax = plt.subplots()

    def get_common_snaps(self):
        ixs_ts_per_sim = [sanitize_tpoints(sd.tpoints) for sd in self.sim_dats]
        # shortest simulation, time wise
        shortest_ix = np.argmin([ix_ts[-1][1] for ix_ts in ixs_ts_per_sim])
        # short simulation final time point
        short_final = ixs_ts_per_sim[shortest_ix][-1][1]
        # ts cropped so that final is <= short_final
        cropped_ixs_ts_per_sim = [
            copy.deepcopy([(ix, t) for (ix, t) in ixs_ts
                           if t < short_final or abs(t - short_final) < 1e-4])
            for ixs_ts in ixs_ts_per_sim
        ]
        # common time points shared by all simulations
        common_ts, snap_ixs_per_sim = find_common_ts(cropped_ixs_ts_per_sim)
        return common_ts, snap_ixs_per_sim

    def investigate(self, d_vertex_plot, d_data_plot, d_max_plot, 
                    d_cell_plot, d_data_group_plot):
        ax.cla()
    
        if abs(d_data_group_plot) > 0:
            DATA_GROUP_IX = (DATA_GROUP_IX + d_data_group_plot) % len(
                DATA_GROUPS)
            ACTIVE_DG = DATA_GROUPS[DATA_GROUP_IX]
            print(ACTIVE_DG.labels)
            NUM_LABEL_GROUPS = len(ACTIVE_DG.label_groups)
            CURR_INNER_IX = 0
    
        CURR_INNER_IX = (CURR_INNER_IX + d_data_plot) % NUM_LABEL_GROUPS
        label_group = ACTIVE_DG.label_groups[CURR_INNER_IX]
    
        VERT_PLOT_IX = (VERT_PLOT_IX + d_vertex_plot) % len(
            VERTEX_PLOT_TYPE)
        CELL_PLOT_IX = (CELL_PLOT_IX + d_cell_plot) % len(CELL_PLOT_TYPE)
    
        vert = VERTEX_PLOT_TYPE[VERT_PLOT_IX]
        cell = CELL_PLOT_TYPE[CELL_PLOT_IX]
    
        if abs(d_max_plot) > 0:
            PLOT_X_MAX = (PLOT_X_MAX + d_max_plot) % (
                    NUM_TSTEPS * NUM_INT_STEPS)
            for g in DATA_GROUPS:
                g.recalc_ylims(PLOT_X_MAX)
    
        for m in range(NUM_CELLS):
            if m == cell or cell == "all":
                py_cell_data = ACTIVE_DG.py_dat[m]
                rust_cell_data = ACTIVE_DG.rust_dat[m]
                if type(label_group) != tuple:
                    tupleized_label_group = (label_group,)
                else:
                    tupleized_label_group = label_group
    
                for label in tupleized_label_group:
                    color = gh.LABEL_COLOR_DICT[label]
                    if len(py_cell_data[label].shape) == 1:
                        ax.plot(
                            rust_cell_data[label][:PLOT_X_MAX],
                            color=color, label=label)
                        ax.plot(
                            py_cell_data[label][:PLOT_X_MAX],
                            color=color,
                            linestyle="dashed", label=label)
                    # ax.set_ylim(ACTIVE_DG.grouped_ylims_dict[label_group])
                    else:
                        for n in range(16):
                            if n == vert or vert == "all":
                                ax.plot(
                                    rust_cell_data[label][:PLOT_X_MAX, n],
                                    color=color, label=label)
                                ax.plot(
                                    py_cell_data[label][:PLOT_X_MAX, n],
                                    color=color,
                                    linestyle="dashed", label=label)
                                # ax.set_ylim(ACTIVE_DG.grouped_ylims_dict[
                                # label_group])
    
        inter_tick_len = np.max([1, np.ceil(PLOT_X_MAX / 20)])
        xticks = np.arange(0, PLOT_X_MAX, inter_tick_len)[:20]
    
        ax.set_xticks(xticks)
        ax.grid(which="major", axis="x")
        ax.legend(loc="best")
        ax.set_title("{}\n{}, vert: {}, cell: {}".format(ACTIVE_DG.description,
                                                         label_group,
                                                         vert, cell))
        fig.canvas.draw()


    def on_press(event):
        global VERT_PLOT_IX
        global CURR_INNER_IX
        global PLOT_X_MAX
        global DATA_GROUP_IX
        global CELL_PLOT_IX
        global fig
        print("pressed: {}".format(event.key))
        if event.key == "down":  # vertex plot change
            paint(-1, 0, 0, 0, 0)
            print(VERT_PLOT_IX)
        elif event.key == "up":  # vertex plot change
            paint(1, 0, 0, 0, 0)
            print(VERT_PLOT_IX)
        elif event.key == "left":  # inner data plot
            paint(0, -1, 0, 0, 0)
            print(CURR_INNER_IX)
        elif event.key == "right":  # inner data plot
            paint(0, 1, 0, 0, 0)
            print(CURR_INNER_IX)
        elif event.key == "z":  # max plot
            paint(0, 0, -1, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "x":  # max plot
            paint(0, 0, 1, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "c":  # max plot
            paint(0, 0, -10, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "v":  # max plot
            paint(0, 0, 10, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "b":  # max plot
            paint(0, 0, -1000, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "n":  # max plot
            paint(0, 0, 1000, 0, 0)
            print(PLOT_X_MAX)
        elif event.key == "home":  # cell plot
            paint(0, 0, 0, 1, 0)
            print(CELL_PLOT_IX)
        elif event.key == "end":  # cell plot
            paint(0, 0, 0, -1, 0)
            print(CELL_PLOT_IX)
        elif event.key == "pagedown":  # data group
            paint(0, 0, 0, 0, -1)
            print(DATA_GROUP_IX)
        elif event.key == "pageup":  # data group
            paint(0, 0, 0, 0, 1)
            print(DATA_GROUP_IX)
        elif event.key == "r":
            VERT_PLOT_IX = 0
            CURR_INNER_IX = 0
            # PLOT_X_MIN = 0
            PLOT_X_MAX = NUM_TSTEPS * NUM_INT_STEPS
            CELL_PLOT_IX = 0
            DATA_GROUP_IX = 0
            paint(0, 0, 0, 0, 0)

    def combined_paint_animation(self, common_t_ix, ax):
        ax.cla()
        ax.set_aspect("equal")
        print("making frame: {}".format(common_t_ix))
        for (sim_ix, sim_dat) in enumerate(self.sim_dats):
            sim_dat.paint_cells(self.snap_ixs_per_sim[sim_ix][common_t_ix], ax)
        ax.set_title("t = {}s".format(self.common_ts[common_t_ix]))
        return ax.get_children()

    def combined_set_ani_opts(self, ani_opts):
        for pls, sim_dat in zip(self.poly_line_styles, self.sim_dats):
            sim_dat.ani_opts = copy.deepcopy(ani_opts)
            sim_dat.ani_opts.poly_line_style = pls

    def animate(self, vec_ani_opts):
        print("beginning combined animation...")
        print("num frames: {}".format(len(self.common_ts)))
        self.ax.set_aspect('equal')
        self.ax.set_xlim(self.default_xlim)
        self.ax.set_ylim(self.default_ylim)
        # shape should be (num_sims, num_snaps)
        for ani_opts in vec_ani_opts:
            self.fig, self.ax = plt.subplots()
            self.combined_set_ani_opts(ani_opts)
            writer = animation.writers['ffmpeg'](fps=10,
                                                 metadata=dict(artist='Me'),
                                                 bitrate=1800)
            cell_ani = animation.FuncAnimation(self.fig,
                                               self.combined_paint_animation,
                                               frames=np.arange(
                                                   len(self.common_ts)),
                                               fargs=(self.ax,),
                                               interval=1, blit=True)
            ani_file_path = self.mp4_file_name + ani_opts.description() + ".mp4"
            ani_save_path = os.path.join(self.out_dir, ani_file_path)
            cell_ani.save(ani_save_path, writer=writer)
            plt.close()
