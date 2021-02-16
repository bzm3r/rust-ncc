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
    rac_act_arrows_per_c_per_s = None
    rho_acts_per_c_per_s = None
    rho_act_arrows_per_c_per_s = None
    adhs_per_c_per_s = None
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

    def load_data(self, out_dir, file_name):
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
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = \
            cb.extract_scalars_from_data(['core', 'rho_acts'], data)
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.adhs_per_c_per_s = \
            cb.extract_p2ds_from_data(['interactions', 'x_adhs'], data)

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
        self.adhs_per_c_per_s = \
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
        self.adhs_per_c_per_s = np.zeros_like(self.poly_per_c_per_s)

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
                self.adhs_per_c_per_s[snap_ix]):
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
