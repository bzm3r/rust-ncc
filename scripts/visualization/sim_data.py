import numpy as np
import matplotlib.pyplot as plt
import os
import cbor2
from matplotlib import animation
from paint_opts import *


def p2ds_to_numpy(p2ds):
    vs = []
    for p2d in p2ds:
        vs.append([p2d['x'], p2d['y']])
    return np.array(vs)


def extract_p2ds_from_cell_states(address, state_recs):
    dat_per_c_per_s = []
    for recs_per_c in state_recs:
        dat_per_c = []
        for cell_rec in recs_per_c:
            relevant_data = cell_rec
            for key in address:
                relevant_data = relevant_data[key]
            dat_per_c.append(p2ds_to_numpy(relevant_data))
        dat_per_c_per_s.append(np.array(dat_per_c))
    return np.array(dat_per_c_per_s)


def extract_scalars(address, state_recs):
    dat_per_c_per_s = []
    for recs_per_c in state_recs:
        dat_per_c = []
        for cell_rec in recs_per_c:
            relevant_data = cell_rec
            for key in address:
                relevant_data = relevant_data[key]
            dat_per_c.append(np.array(relevant_data))
        dat_per_c_per_s.append(np.array(dat_per_c))
    return np.array(dat_per_c_per_s)


class SimulationData:
    def __init__(self, out_dir, file_name):
        self.out_dir = out_dir
        self.cbor_file_path = file_name + ".cbor"
        self.mp4_file_name = file_name

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
                "Error: could not find requested file {} in dir with "
                "contents: {}".format(self.cbor_file_path, cbor_files))
        self.out_file_path = os.path.join(self.out_dir, self.cbor_file_path)

        snapshots = []
        with open(self.out_file_path, mode='rb') as sf:
            world_info = cbor2.load(sf)
            success = True
            while success:
                try:
                    snapshots += cbor2.load(sf)
                finally:
                    success = False

        print("file_name: {} | snapshots found: {}".format(file_name,
                                                           len(snapshots)))
        self.tsteps = [s["tstep"] for s in snapshots]
        state_recs = [s["cells"] for s in snapshots]
        self.frequency = world_info["snap_period"]

        self.poly_per_c_per_s = \
            extract_p2ds_from_cell_states(['core', 'poly'], state_recs)
        self.centroids_per_c_per_s = np.array(
            [[np.average(poly, axis=0) for poly in poly_per_c] for
             poly_per_c in
             self.poly_per_c_per_s])
        self.uivs_per_c_per_s = \
            extract_p2ds_from_cell_states(['core', 'geom', 'unit_in_vecs'],
                                          state_recs)
        self.uovs_per_c_per_s = -1 * self.uivs_per_c_per_s
        self.rac_acts_per_c_per_s = extract_scalars(['core', 'rac_acts'],
                                                    state_recs)
        self.rac_act_arrows_per_c_per_s = \
            self.rac_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uovs_per_c_per_s
        self.rho_acts_per_c_per_s = extract_scalars(['core', 'rho_acts'],
                                                    state_recs)
        self.rho_act_arrows_per_c_per_s = \
            self.rho_acts_per_c_per_s[:, :, :, np.newaxis] * \
            self.uivs_per_c_per_s
        self.adhs_per_c_per_s = \
            extract_p2ds_from_cell_states(['interactions', 'x_adhs'],
                                          state_recs)
        self.centroids_per_c_per_s = np.array(
            [[np.average(poly, axis=0) for poly in poly_per_cell] for
             poly_per_cell in
             self.poly_per_c_per_s])

        self.snap_ix = 0
        self.default_xlim = None
        self.default_ylim = None
        self.default_bbox_lim = None
        self.ani_opts = default_ani_opts()
        self.fig_probe, self.ax_probe = None, None
        self.fig_ani, ax = None, None

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
        self.fig_ani, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(self.default_xlim)
        ax.set_ylim(self.default_ylim)
        frame_ixs = [n for n in range(int(len(self.tsteps)))]
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
                                               fargs=(ax,),
                                               interval=1, blit=True)
            ani_file_path = self.mp4_file_name + ani_opts.description() + ".mp4"
            ani_save_path = os.path.join(self.out_dir, ani_file_path)
            cell_ani.save(ani_save_path, writer=writer)

    def paint_cells(self, snap_ix, ax):
        for (ci, poly) in enumerate(self.poly_per_c_per_s[snap_ix]):
            for vix in range(16):
                ax.plot([poly[vix, 0], poly[(vix + 1) % 16, 0]],
                        [poly[vix, 1], poly[(vix + 1) % 16, 1]],
                        color="k", marker=".", markersize="0.5")
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
                self.tsteps[self.snap_ix], self.snap_ix
            )
        )
        self.snap_ix = (self.snap_ix + delta) % len(self.tsteps)
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
        ax.relim()

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

        ax.set_title(
            "tstep {} (snapshot {})".format(
                self.tsteps[snap_ix], snap_ix
            )
        )
        return ax.get_children()
