import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from paint_opts import *
import get_comp as cd
import get_cbor as cb
from utils import *
from sim_data import *


class PythonRustComparisonData:
    def __init__(self, out_dir, py_dat, rust_dat, poly_line_styles,
                 mp4_file_name):
        self.out_dir = out_dir
        self.py_dat = py_dat
        self.rust_dat = rust_dat
        self.sim_dats = [py_dat, rust_dat]
        self.verify_parameter_equality()
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

    def verify_parameter_equality(self):
        rust_keys = self.rust_dat.header.keys()
        py_keys = self.py_dat.header.keys()

        not_in_py = []
        for key in rust_keys:
            if key not in py_keys:
                not_in_py.append(key)

        not_in_rust = []
        for key in py_keys:
            if key not in rust_keys:
                not_in_rust.append(key)

        print("not_in_py: {}".format(not_in_py))
        print("not_in_rust: {}".format(not_in_rust))

        for key in rust_keys:
            if key in py_keys:
                rust_param = self.rust_dat.header[key]
                py_param = self.py_dat.header[key]
                print("{}: rust = {}, py = {}".format(key, rust_param,
                                                      py_param))
                if type(rust_param) == list:
                    rust_param = np.array(rust_param)
                    py_param = np.array(py_param)
                    delta = np.max(abs(rust_param - py_param))
                else:
                    delta = abs(rust_param - py_param)
                if delta > 1e-4:
                    raise Exception(
                        "parameter mismatch {}: rust = {}, py = {}. delta = {}"
                            .format(key, self.rust_dat.header[key],
                                    self.py_dat.header[key], delta))

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

    def combined_set_ani_opts(self, ani_opts):
        for pls, sim_dat in zip(self.poly_line_styles, self.sim_dats):
            sim_dat.ani_opts = copy.deepcopy(ani_opts)
            sim_dat.ani_opts.poly_line_style = pls

    def combined_paint_animation(self, common_t_ix, ax, ty):
        ax.cla()
        ax.set_aspect("equal")
        print("making frame: {}".format(common_t_ix))
        for (sim_ix, sim_dat) in enumerate(self.sim_dats):
            sim_dat.paint_cells(self.snap_ixs_per_sim[sim_ix][common_t_ix],
                                ax, ty)
        ax.set_title("t = {}s".format(self.common_ts[common_t_ix]))
        return ax.get_children()

    def animate(self, vec_ani_opts, ty):
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
                                               fargs=(self.ax, ty),
                                               interval=1, blit=True)

            if len(ty) != 0:
                ty_tag = "_{}".format(ty)
            else:
                ty_tag = ""
            ani_file_path = self.mp4_file_name + ani_opts.description() + \
                            ty_tag + ".mp4"
            ani_save_path = os.path.join(self.out_dir, ani_file_path)
            cell_ani.save(ani_save_path, writer=writer)
            plt.close()

    def plot(self, data_labels):
        rs = [eval("self.rust_dat.{}_per_c_per_s".format(dl), {"self": self})
              for
              dl in
              data_labels]
        ps = [eval("self.py_dat.{}_per_c_per_s".format(dl), {"self": self}) for dl in
              data_labels]

        for r, p, l in zip(rs, ps, data_labels):
            fig, ax = plt.subplots()
            cell_ix = 0
            ax.plot(p[:, cell_ix, :1], label="python", linewidth=3, alpha=0.3)
            ax.plot(r[:, cell_ix, :1], label="rust", linewidth=1, alpha=1.0)
            ax.set_title("{} for cell {}".format(l, cell_ix))
            ax.legend(loc="best")
            plot_path = os.path.join(self.out_dir, "{}_{}.png".format(
                self.mp4_file_name, l))
            fig.savefig(plot_path)
            plt.close()
