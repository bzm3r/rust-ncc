from retrieve import *
import numpy as np
import matplotlib.pyplot as plt
import copy
import py_model.hardio as hio
import graph_help as gh

NUM_TSTEPS = 18
NUM_CELLS = 2
COA = 0
CIL = 0
NUM_INT_STEPS = 10

py_out = read_save_file("py", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL, COA)
rust_out = read_save_file("rust", NUM_TSTEPS, NUM_INT_STEPS, NUM_CELLS, CIL,
                          COA)

check_header_equality(py_out, rust_out)

py_dat = py_out["tsteps"]
rust_dat = rust_out["tsteps"]

py_data_dict_per_cell = gen_data_dict_per_cell_int_steps(py_out)
rust_data_dict_per_cell = gen_data_dict_per_cell_int_steps(rust_out)


def paint(delta_vertex_plot, delta_data_plot, delta_max_plot,
          delta_cell_plot, delta_data_group_plot):
    global fig
    global ax
    global VERT_PLOT_IX
    global CURR_INNER_IX
    global PLOT_X_MAX
    global DATA_GROUP_IX
    global CELL_PLOT_IX
    global NUM_CELLS
    global DATA_GROUPS
    global ACTIVE_DG
    global NUM_LABEL_GROUPS
    global LABEL_GROUPS
    global DG_COLORS
    global DG_YLIMS
    global GROUP_PY_DATA
    global GROUP_RUST_DATA
    global GROUP_DESCRIPT
    ax.cla()

    if abs(delta_data_group_plot) > 0:
        DATA_GROUP_IX = (DATA_GROUP_IX + delta_data_group_plot) % len(
            DATA_GROUPS)
        CURR_DATA_GROUP = DATA_GROUPS[DATA_GROUP_IX]
        GROUP_DATA_IX_MAX = len(CURR_DATA_GROUP.labels)
        GROUP_LABELS = CURR_DATA_GROUP.labels
        GROUP_COLORS = [gh.LABEL_COLOR_DICT[label] for label in GROUP_LABELS]
        GROUP_YLIMS = CURR_DATA_GROUP.grouped_ylims_dict
        GROUP_PY_DATA = CURR_DATA_GROUP.py_dat
        GROUP_RUST_DATA = CURR_DATA_GROUP.rust_dat
        GROUP_DESCRIPT = DATA_GROUPS[DATA_GROUP_IX].description
        CURR_DATA_IX = 0

    CURR_DATA_IX = (CURR_DATA_IX + delta_data_plot) % GROUP_DATA_IX_MAX
    label = GROUP_LABELS[CURR_DATA_IX]
    color = GROUP_COLORS[CURR_DATA_IX]

    VERTEX_IX = (VERTEX_IX + delta_vertex_plot) % len(
        VERTEX_PLOT_TYPE)
    CELL_IX = (CELL_IX + delta_cell_plot) % len(CELL_PLOT_TYPE)

    vert = VERTEX_PLOT_TYPE[VERTEX_IX]
    cell = CELL_PLOT_TYPE[CELL_IX]

    if abs(delta_max_plot) > 0:
        PLOT_X_MAX = (PLOT_X_MAX + delta_max_plot) % (NUM_TSTEPS * NUM_INT_STEPS)
        for g in DATA_GROUPS:
            g.recalc_ylims(PLOT_X_MAX)
        GROUP_YLIMS = CURR_DATA_GROUP.grouped_ylims_dict
        GROUP_PY_DATA = CURR_DATA_GROUP.py_dat
        GROUP_RUST_DATA = CURR_DATA_GROUP.rust_dat

    for m in range(NUM_CELLS):
        if m == cell or cell == "all":
            dict_py_dat = GROUP_PY_DATA[m]
            dict_rust_dat = GROUP_RUST_DATA[m]
            if len(dict_rust_dat[label].shape) == 1:
                ax.plot(
                    dict_rust_dat[label][:PLOT_X_MAX] - dict_py_dat[label][
                                                        :PLOT_X_MAX],
                    color=color)
                # ax.plot(
                #     dict_py_dat[label][:PLOT_X_MAX],
                #     color=color,
                #     linestyle="dashed")
                # ax.set_ylim(GROUP_YLIMS[label])
            else:
                for n in range(16):
                    if n == vert or vert == "all":
                        ax.plot(dict_rust_dat[label][
                                :PLOT_X_MAX, n] - dict_py_dat[label][
                                                     :PLOT_X_MAX, n],
                                color=color)
                        # ax.plot(dict_py_dat[label][
                        #         :PLOT_X_MAX, n], color=color,
                        #         linestyle="dashed")
                        print(label)
                        print(GROUP_YLIMS)
                        ylims = GROUP_YLIMS[label]
                        # ax.set_ylim(ylims)

    # ax.set_xticks(np.arange(PLOT_X_MAX))
    # ax.grid(which="major", axis="x")
    ax.set_title("{}\n{}, vert: {}, cell: {}".format(GROUP_DESCRIPT, label,
                                                     vert, cell))


def on_press(event):
    global VERT_PLOT_IX
    global CURR_INNER_IX
    global PLOT_X_MAX
    global DATA_GROUP_IX
    global CELL_PLOT_IX
    global fig
    print("pressed: {}".format(event.key))
    if event.key == "up":
        paint(1, 0, 0, 0, 0)
    elif event.key == "down":
        paint(-1, 0, 0, 0, 0)
    elif event.key == "right":
        paint(0, 1, 0, 0, 0)
    elif event.key == "left":
        paint(0, -1, 0, 0, 0)
    elif event.key == "z":
        paint(0, 0, -10, 0, 0)
    elif event.key == "x":
        paint(0, 0, 10, 0, 0)
    elif event.key == "c":
        paint(0, 0, -1, 0, 0)
    elif event.key == "v":
        paint(0, 0, 1, 0, 0)
    elif event.key == " ":
        paint(0, 0, 0, 1, 0)
    elif event.key == "pagedown":
        paint(0, 0, 0, 0, -1)
    elif event.key == "pageup":
        paint(0, 0, 0, 0, 1)
    elif event.key == "r":
        VERTEX_IX = 0
        CURR_DATA_IX = 0
        PLOT_X_MAX = 150
        CELL_IX = 0
        DATA_GROUP_IX = 0
        paint(0, 0, 0, 0, 0)


PLOT_X_MAX = int(0.01 * NUM_TSTEPS * NUM_INT_STEPS)
VERTEX_PLOT_TYPE = [n for n in range(16)] + ["all"]
CELL_PLOT_TYPE = [m for m in range(NUM_CELLS)] + ["all"]
VERT_PLOT_IX = 0
CELL_PLOT_IX = 0

py_data_dict_per_cell, rust_data_dict_per_cell, DATA_GROUPS = \
    gh.generate_groups(
        py_data_dict_per_cell, rust_data_dict_per_cell, PLOT_X_MAX)
DATA_GROUP_IX = 0
ACTIVE_DG = DATA_GROUPS[DATA_GROUP_IX]
NUM_LABEL_GROUPS = len(ACTIVE_DG.labels)
LABEL_GROUPS = ACTIVE_DG.labels
DG_COLORS = [gh.LABEL_COLOR_DICT[label] for label in LABEL_GROUPS]
DG_YLIMS = ACTIVE_DG.grouped_ylims_dict
GROUP_PY_DATA = ACTIVE_DG.py_dat
GROUP_RUST_DATA = ACTIVE_DG.rust_dat
GROUP_DESCRIPT = DATA_GROUPS[DATA_GROUP_IX].description
CURR_INNER_IX = 0

fig, ax = plt.subplots()
paint(0, 0, 0, 0, 0)
fig.canvas.mpl_connect('key_press_event', on_press)
